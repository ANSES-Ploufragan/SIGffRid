#!/usr/bin/perl
use strict;
# use lib '/Volumes/applications/Applications/bin/bioperl-1.5.2_102/';
use Env qw(SIGFFRID_DIR);
use Bio::SeqIO;
use Getopt::Long;
use Cwd;

# ************************************************************************************
# DECLARE STORING VARIABLES AND PARAMETERS
# ************************************************************************************

use warnings;
#use diagnostics;

# personal LIB ********************************************************************************
use lib $SIGFFRID_DIR.'libperl';
use SysUtilities; # qw(date tri_rapide)
use SeqUtilities; # qw(complement_seq print50 get_cosmid_seq);
# END LIB *************************************************************************************

# **********************************************************************

=head1  NAME

FextractSeqFeaturesFromGbUsingBioPerl.pl

=head1 DESCRIPTION

Extract sequences from gb or gbk files with various properties.
Mainly dedicated to dataset for SIGffRid analysis.
Caution, if output file already exists, overwrite it.

=head1  USAGE

=over

=item -gb_f <s>n

Genbank file for bacteria (must contain nucleotide sequence).

=item -id_bact <s>

ID of bacteria (this id is in the file provided by MBGD database
file of ortholog comparisons.

=item [-res_dir <s>]

Default:current directory.

=item [-test]

To run only test of seq extraction, do not take into account other arguments.

=back

examples:


Will run for example following command lines:

perl sigffrid_cmd_line.pl -gb_f \
../ref_genomes/NC_002516.2_Pseudomonas_aeruginosa_PAO1_chr_complete_genome.gb
-id_bact pae
-res_dir sigffrid_res/


=cut

# **********************************************************************

my $prog_tag = '[FextractSeqFeaturesFromGbUsingBioPerl]';
my $b_test   = 0; # ok 2018 12 14
my @gb_f     = ();
my $res_dir  = undef;

# we get date with '/' characte to separate day, month and year 
my $date = &date('/');
print "$prog_tag DATE: $date\n";

# booleens de choix ou non d'un type de traitement (dans l'ordre d'apparation dans le traitement)
my $debug_mode = 0;

my $ID_bact  = undef;

# format GERE DANS LA FONCTION PRINT50
my $fasta                 = 0; # edition des sequences au format fasta 1 ou comme une seul sequence sur une ligne pour accelerer le traitement de certains programmes.

my $entete_fasta;              # en tete de fichier pour preciser le format utilise pour generer le fichier
if($fasta){ $entete_fasta = 'fasta_'     }
else      { $entete_fasta = 'not_fasta_' }

my $edit_file_for_the_whole_one_block_seq0      = 1; # edition de la sequence complete d'une bacterie en sens 0:  OK
my $edit_file_for_the_whole_one_block_seq1      = 0; # edition de la sequence complete d'une bacterie en sens 1:  OK
my $edit_file_for_the_whole_one_block_seq_RMES  = 1; # edition de la sequence complete d'une bacterie en sens 0, si plusieurs seq, separees par Z, une entete obligatoire
my $edit_file_for_each_CDS                      = 0; # edition de chaque CDS individuellement (meme en cas de chevauchement):  OK
my $edit_file_for_all_CDS                       = 0; # edition de l'ensemble des CDS (fusions en cas de chevauchement): OK
my $edit_file_for_promoters                     = 1; # edition du promoteur de chaque gene (seq interg amont du CDS): OK
my $edit_file_for_all_interCDS                  = 0; # edition de l'ensemble des interCDS
my $edit_file_for_all_interCDS_between_conv     = 0; # edition de l'ensemble des interCDS entre genes conv
my $edit_file_for_all_interCDS_between_div      = 0; # edition de l'ensemble des interCDS entre genes div
my $edit_file_for_all_CDS_interCDS_between_conv = 0; # edition de l'ensemble des CDS+interCDS entre genes conv : prev OK + limit
my $edit_file_for_all_annotations               = 0; # edition de toutes les annotations correspondant a chaque CDS OK
my $edit_file_for_upstream_seq                  = 1; # edition de la seq amont de chaque gene (seq amont du CDS), seq CDS et/ou interCDS : OK
my $edit_file_for_protein_seq                   = 0; # edition de toutes les sequences proteiques: OK 
my $edit_file_for_300_0_seq                     = 1; # 1; # edition de toutes les sequences -300 0 en sens 0 puis en sens 1 (fusion pour celle se chevauchant dans le meme sens, sans chevauchement pour un sens considere: OK
my $edit_limits_of_overlapped_seq_sens01        = 1; # 1; # edition de toutes les limites des séquences amonts fusionnées: OK
my $edit_annotations_of_limits_overlapped       = 1; # 1; # edition des annotations des genes dont les limites amonts sont fusionnées: OK
my $edit_limits_of_promoters                    = 0; # future work... edition des limites des sequences intergéniques amonts:OK

# limites utilisees pour l'extraction des sequences amont de CDS (qu'elles soient codantes ou non)
my $longueurmini             = 30; # longueur minimale de la sequence extraite
my $dist_fin_up_seq_deb_trad = 0; # distance de la fin de la sequence amont au debut de la traduction
my $dist_deb_up_seq_deb_trad = 350; # distance a partir de laquelle on considere la sequence amont
my $length_upstream = $dist_deb_up_seq_deb_trad - $dist_fin_up_seq_deb_trad; # position du debut de la seq amont - position de fin

my $not_verbose                                 = 1; # affichage des infos sur le deroulement du prog si = 0
my $bool_extended_alphabet                      = 1; # utilisation de l'alphabet etendu pour les seq qui ont d'autres lettres que acgt

my $bool_distinctfiles_sens0_1 = 0; # to delete files of upstream seq sens 0 and sens 1 after creatinf file sens01
# VAR for write in output files ***************************************
my @tmp_file_name          = (); # store output files names
my $last_ind_1sens_file    = -1;

my $end_file = '';

my $chemin_dest = undef;

# **********************************************************************
# CHECK OPTIONS
# **********************************************************************
my $nbargmini = 1;
if(scalar(@ARGV) < $nbargmini){ 
  print "$prog_tag Bad number of arguments: ".scalar(@ARGV)." found, at least $nbargmini wanted\n";
  foreach(0..$#ARGV){
    print "$prog_tag $_: $ARGV[$_]\n";
  }
  die `perldoc $0`;
}
GetOptions(
    "gb_f=s{1,}"              => \@gb_f,
    "id_bact=s"               => \$ID_bact,
    "res_dir=s"               => \$chemin_dest,
    "test"                    => sub { $b_test  = 1 }
);
# **********************************************************************

if($b_test)
{
    $ID_bact = 'H37Rv';
    my $test_dir = $SIGFFRID_DIR.'test_SIGffRid/test_FextractSeqFeaturesFromGbUsingBioPerl/';
    @gb_f = ($test_dir.'AL123456.3_Mycobacterium_tuberculosis_H37Rv_complete_genome_full.gb');    
    $res_dir = $test_dir.'res/';
    my $cmd = join(' ', $SIGFFRID_DIR.'FextractSeqFeaturesFromGbUsingBioPerl.pl',
		   '-gb_f', @gb_f,
		   '-id_bact', $ID_bact,
		   '-res_dir', $res_dir
	);
    print "$prog_tag [TEST] cmd:$cmd, line ".__LINE__."\n";
    `$cmd`;
    print "$prog_tag [TEST] end\n";
    exit;
}
# if($#ARGV < 1){ &sortir;}
if( defined $chemin_dest )
{  
    if($chemin_dest !~ /\/$/)
    {	$chemin_dest .= '/';    }
    -e $chemin_dest or die "$prog_tag [Error] $chemin_dest does not exist, line ".__LINE__."\n";
    -d $chemin_dest or die "$prog_tag [Error] $chemin_dest is not a directory, line ".__LINE__."\n";    
}
else
{    
  $chemin_dest = Cwd::cwd(); 
  $chemin_dest .= '/';
}
foreach(@gb_f){
  my $b_sortir = 0;
  if(!-e $_){
    die "$prog_tag [Error] $_ file does not exist, line ".__LINE__."\n";
  }
  if($_ !~ /\.gbk?$/){
    die "$prog_tag [Error] $_ does not have GenBank file extension (.gb or .gbk), line ".__LINE__."\n";
  }
}

# sub sortir{
#  print "Usage:$0 <ID bacteria> <GenBank file(s) (file.gbk or file.gb)> [destination directory]\n";
#  exit;
# }

print "$prog_tag chemin_dest : $chemin_dest\n";
# my $ID_bact  = $ARGV[0] ;

my %H_ind_T = (
	       'gene'=> 0,
	       'prot'=>1,
	       'imp'=>2
	      );
my $bool_premier_passage = 1;
my $cpt_seq = 0;
my $separateur = '_';

my $verbose       = 0;        # to display or not debug informations
my  @gene          = ();       # store genes names (coding prot OR rna)   
my  @comment       = ();       # store comments for genes (prot and rna), imp, cdregion...
my  @seqdesc_title = ();
my  @sens          = ();       # store strand (plus = 0, minus = 1, both = 2, unknown = 3) for genes or rna [0], prot [1], imp [2]
my  @prot          = ();       # store prot seq
my  @prot_name     = ();       # store prot name
my  @Nr_acc_prot   = ();       # store prot IDs (NP_, ZP_...)
my  @type          = ();       # store type (prot, tRNA...)
my  @beg = ();       # store beginning (limit) of CDS [0], prot [1], imp [2], ... (do not take into account sens of gene)
my  @end = ();       # store end (limit) of CDS (do not take into account sens of gene)
# foreach(keys %H_ind_T){
#   @{ $beg[ $H_ind_T{$_} ] } = ();
#   @{ $end[ $H_ind_T{$_} ] } = ();
#   @{ $sens[ $H_ind_T{$_} ] } = ();
# }

# ************************************************************************************
# END DECLARE STORING VARIABLES AND PARAMETERS
# ************************************************************************************

# ************************************************************************************
# EXTRACTION
# ************************************************************************************

# &GetOptions("r"=>\$reverse);


for my $gbkfile(@gb_f) {

  @gene          = ();       # store genes names (coding prot OR rna)   
  @comment       = ();       # store comments for genes (prot and rna), imp, cdregion...
  @seqdesc_title = ();
  @sens          = ();       # store strand (plus = 0, minus = 1, both = 2, unknown = 3) for genes or rna [0], prot [1], imp [2]
  @prot          = ();       # store prot seq
  @prot_name     = ();       # store prot name
  @Nr_acc_prot   = ();       # store prot IDs (NP_, ZP_...)
  @type          = ();       # store type (prot, tRNA...)
  @beg = ();       # store beginning (limit) of CDS [0], prot [1], imp [2], ... (do not take into account sens of gene)
  @end = ();       # store end (limit) of CDS (do not take into account sens of gene)

  my $sequence = &get_cosmid_seq($gbkfile,$bool_extended_alphabet);
  my $length_sequence = length($sequence);
  print "$prog_tag length_sequence : $length_sequence\n";

  my $k      = 0;
  my $length = 0;

  print "$prog_tag Treatment of $gbkfile file\n";

  my $seqIO_obj=Bio::SeqIO->new(-format=>'GenBank', -file=>"$gbkfile");
  
  my $pb_in_genes=0;

  while(my $seq_obj=$seqIO_obj->next_seq()){
    $length=$seq_obj->length();
    my @features=$seq_obj->all_SeqFeatures();

    foreach my $feature (@features){
      my $pri_tag=$feature->primary_tag();
      unless (($pri_tag eq 'CDS')||($pri_tag =~ /RNA/)) {next;}
      
      my @genes_list = ();
      my @locus_list = ();
      my @product    = ();
      my @note       = ();
      my @protein_id = ();
      my $bool_found = 0;
      my $genes = '';
      if($feature->has_tag('gene')){
	@genes_list=$feature->each_tag_value('gene');
	$bool_found = 1;
      }
      if($feature->has_tag('locus_tag')){
	@locus_list=$feature->each_tag_value('locus_tag');
	$bool_found = 1;
      }
      if($bool_found){

	$genes=join('_',@genes_list,@locus_list);
	$genes=~s/\ //g;
      }
      else{ $genes=$pri_tag; }
      
      if($feature->has_tag('product')){
	@product=$feature->get_tag_values('product');
      }
      
      if($feature->has_tag('note')){
	@note=$feature->get_tag_values('note');
      }
      
      if($feature->has_tag('protein_id')){
	@protein_id=$feature->get_tag_values('protein_id');
      }
      
      my $start=$feature->start();
      my $end=$feature->end();
      my $strand=$feature->strand();
      $k++;
      $pb_in_genes+=$end-$start+1;
      print("$genes\t$start\t$end\t$strand\t@product\t@protein_id\n"); # \t@note\n");
      
      # RECORD IN DATA STRUCT

	push @gene, $genes;
	push @{ $beg[ $H_ind_T{'gene'} ] }, $start;
	push @{ $end[ $H_ind_T{'gene'} ] }, $end;
	push @{ $sens[ $H_ind_T{'gene'} ] }, $strand;
	push @prot_name, join(' ',@product);
	push @Nr_acc_prot, join(' ',@protein_id);
	push @{ $comment[ $H_ind_T{'gene'} ] }, join(' ',@note);
	push @type, $pri_tag;

    }
  }

  print("parsed $k genes in $length pb sequence with $pb_in_genes pb in genes and ",$length-$pb_in_genes," pb in intergenic spaces (");
  printf("%.2f",100*($length-$pb_in_genes)/$length);
  print "\%)\n";

# ************************************************************************************
# END EXTRACTION 
# ************************************************************************************


# ************************************************************************************
# WRITING IN FILES
# ************************************************************************************

   $date = &date();

# record of the whole seq in sens 1 (plus)
if( $edit_file_for_the_whole_one_block_seq0 ){
    
    if($bool_premier_passage and (-e "$chemin_dest".$entete_fasta.'whole_seq0_'.$ID_bact."$end_file".".txt") )
    {
    	system("rm -f $chemin_dest".$entete_fasta.'whole_seq0_'.$ID_bact."$end_file".".txt");
    }

    open(FSORWHOLESEQ0,">> $chemin_dest".$entete_fasta.'whole_seq0_'.$ID_bact."$end_file".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".$entete_fasta."whole_seq0_".$ID_bact."$end_file".".txt\n";
    print FSORWHOLESEQ0 &print50($fasta,$sequence);
    close FSORWHOLESEQ0;
   
    print "$prog_tag length_sequence : ".length($sequence)."\n";
    
    print "$prog_tag $chemin_dest".$entete_fasta.'whole_seq0_'.$ID_bact."$end_file"."_$date".".txt file created\n";
}

# record of the whole seq in sens -1 (minus)
if( $edit_file_for_the_whole_one_block_seq1 ){
    my $rev_seq = reverse "$sequence";
    $rev_seq =~ tr/acgtACGT/tgcaTGCA/;

    if($bool_premier_passage and (-e "$chemin_dest".$entete_fasta.'whole_seq1_'.$ID_bact."$end_file".".txt") )
    {
   	system("rm -f $chemin_dest".$entete_fasta.'whole_seq1_'.$ID_bact."$end_file".".txt");
    }

    open(FSORWHOLESEQ1,">> $chemin_dest".$entete_fasta.'whole_seq1_'.$ID_bact."$end_file".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".$entete_fasta."whole_seq1_".$ID_bact."$end_file".".txt\n";
    # &complement_seq(\$sequence);
    print FSORWHOLESEQ1 &print50($fasta,$rev_seq);
    close FSORWHOLESEQ1;
    
    print "$prog_tag $chemin_dest".$entete_fasta.'whole_seq1_'.$ID_bact."$end_file"."_$date".".txt file created\n";
}

# edit whole seq(s) for rmes prog
if($edit_file_for_the_whole_one_block_seq_RMES){

#     if(($bool_premier_passage) and (-e $chemin_dest.'RMES_'.$ID_bact.$end_file.".txt") )
#     {
# 	die "File $chemin_dest".'RMES_'.$ID_bact.$end_file.".txt already exists, change the destination folder \n";
#     }
#     
    if($bool_premier_passage and (-e "$chemin_dest".'RMES_'.$ID_bact.$end_file.".txt"))
    {
	system("rm -f $chemin_dest".'RMES_'.$ID_bact.$end_file.".txt"); 
    }

    open(FRMES,">> $chemin_dest".'RMES_'.$ID_bact.$end_file.".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".'RMES_'.$ID_bact.$end_file.".txt\n";

    if($bool_premier_passage){ print FRMES ">$ID_bact\n"; }
    else                     {  print FRMES "ZZZZZZZZZZ";        }

    print FRMES $sequence;
    close FRMES;
   
    print "$prog_tag length_sequence : ".length($sequence)."\n";
    
    print $prog_tag.' '.$chemin_dest.'RMES_'.$ID_bact.$end_file.".txt file created\n";
}

# record of every CDS sequences
if($edit_file_for_each_CDS){

#     if(($bool_premier_passage) and (-e "$chemin_dest".$entete_fasta.'chaque_CDS_'.$ID_bact."$end_file"."_$date".".txt") )
#     {
# 	die "File $chemin_dest".$entete_fasta.'chaque_CDS_'.$ID_bact."$end_file"."_$date".".txt already exists, change the destination folder \n";
#     }

    if($bool_premier_passage and (-e "$chemin_dest".$entete_fasta.'chaque_CDS_'.$ID_bact."$end_file"."_$date".".txt"))
    {
	system("rm -f $chemin_dest".$entete_fasta.'chaque_CDS_'.$ID_bact."$end_file"."_$date".".txt");    
    }
    # pour toutes les sequences CDS (chevauchements possibles)
    open(FSOREACHCDS,">> $chemin_dest".$entete_fasta.'chaque_CDS_'.$ID_bact."$end_file"."_$date".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".$entete_fasta."chaque_CDS_".$ID_bact."$end_file"."_$date".".txt\n";


    for my $cpt_deb_CDS(0..$#{$beg[ $H_ind_T{'gene'} ]})
    {

	    # if($beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS]+1 < $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS])
	    # { 
	    # enregistrement du codant correspondant
	    my $long = $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] +1;
	    print FSOREACHCDS ">CDS_sequence_of_$gene[$cpt_deb_CDS], sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] limits : $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS],$end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]\n";    
	    ($not_verbose)or print "\n>CDS sequence of $gene[$cpt_deb_CDS] sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS]\n";
	    if($sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == -1)
	    {
		print FSOREACHCDS &print50($fasta,&complement_seq(substr($sequence,$beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] -1, $long)));
		($not_verbose)or print &print50($fasta,&complement_seq(substr($sequence,$beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] -1, $long)));
	    }
	    elsif($sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == 1)
	    {
		print FSOREACHCDS &print50($fasta,substr($sequence,$beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] -1, $long));
		($not_verbose)or print &print50($fasta,substr($sequence,$beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] -1, $long));
	    }
	    else
	    {
		print "$prog_tag Unknown sens for $gene[$cpt_deb_CDS]\n";
	    }
	    # }
	    # else
	    # {
	    # die "Problem : end of CDS ($end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]) lower than beginning ($beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS]) for gene $gene[$cpt_deb_CDS]\n";
	    # }

    }
    close FSOREACHCDS;
    print "$prog_tag $chemin_dest".$entete_fasta.'chaque_CDS_'.$ID_bact."$end_file"."_$date".".txt file created\n";
}

if($edit_file_for_all_CDS)
{
#    if(($bool_premier_passage) and (-e "$chemin_dest".$entete_fasta.'ens_seq_CDS_'.$ID_bact."$end_file"."_$date".".txt") )
#     {
# 	die "File $chemin_dest".$entete_fasta.'ens_seq_CDS_'.$ID_bact."$end_file"."_$date".".txt already exists, change the destination folder \n";
#     }
#     

    if($bool_premier_passage and (-e "$chemin_dest".$entete_fasta.'ens_seq_CDS_'.$ID_bact."$end_file"."_$date".".txt")){
	system("rm -f $chemin_dest".$entete_fasta.'ens_seq_CDS_'.$ID_bact."$end_file"."_$date".".txt");
    }
# pour l'ensemble des sequences CDS (pas de chevauchement)
    open(FSORENSCDS,">> $chemin_dest".$entete_fasta.'ens_seq_CDS_'.$ID_bact."$end_file"."_$date".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".$entete_fasta."ens_seq_CDS_".$ID_bact."$end_file"."_$date".".txt\n";
    
    ($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
    my ($debut, $fin) = &get_limits(0);
    my ($gene_beg, $gene_end);
    my $temp_deb  = undef;
    for my $cpt_deb_CDS(1..$#{$beg[ $H_ind_T{'gene'} ]})
    {
	($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
	($gene_beg, $gene_end) = &get_limits($cpt_deb_CDS);
	if($gene_beg > $fin+1)
	{ 
# enregistrement du codant precedent
	    $cpt_seq++;
	    my $long = $fin - $debut +1;
	    print FSORENSCDS ">CDS_sequence_set_$cpt_seq limits : $debut, $fin\n";
	    print FSORENSCDS &print50($fasta,substr($sequence,$debut -1, $long));
	    ($not_verbose)or print "\n>CDS Sequence set $cpt_seq\n";
	    ($not_verbose)or print &print50($fasta,substr($sequence,$debut -1, $long));
	    ($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
	    ($debut, $fin) = &get_limits($cpt_deb_CDS);
	    if($cpt_deb_CDS == $#{$beg[ $H_ind_T{'gene'} ]})
	    {
		print FSORENSCDS ">CDS_sequence_set_$cpt_seq limits : $debut, $fin\n";
		print FSORENSCDS &print50($fasta,substr($sequence,$debut -1, $long));
	    }
	   
	}
	else
	{
	    ($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
	    ($temp_deb, $fin) = &get_limits($cpt_deb_CDS);
	}
    }
    close FSORENSCDS;

    print "$prog_tag $chemin_dest".$entete_fasta.'ens_seq_CDS_'.$ID_bact."$end_file"."_$date".".txt file created\n";
    $cpt_seq = 0;
}


# ******** ??? A REPRENDRE
if($edit_file_for_promoters)
{
#     if(($bool_premier_passage) and (-e "$chemin_dest".$entete_fasta.'promot_'.$ID_bact."$end_file".".txt") )
#     {
# 	die "File $chemin_dest".$entete_fasta.'promot_'.$ID_bact."$end_file".".txt already exists, change the destination folder \n";
#     }

    if($bool_premier_passage and (-e "$chemin_dest".$entete_fasta.'promot_'.$ID_bact."$end_file".".txt")){ 
	system("rm -f $chemin_dest".$entete_fasta.'promot_'.$ID_bact."$end_file".".txt"); 
    }
# output file pour toutes les sequences promotrices (chevauchements possibles)
    open(FSORPROMOT,">> $chemin_dest".$entete_fasta.'promot_'.$ID_bact."$end_file".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".$entete_fasta."promot_".$ID_bact."$end_file".".txt\n";
    if($edit_limits_of_promoters)
    {
# 	if(($bool_premier_passage) and (-e "$chemin_dest".'limits_of_promot_'.$ID_bact."$end_file".".txt") )
# 	{
# 	    die "File $chemin_dest".'limits_of_promot_'.$ID_bact."$end_file".".txt already exists, change the destination folder \n";
# 	}
	if($bool_premier_passage and (-e "$chemin_dest".'limits_of_promot_'.$ID_bact."$end_file".".txt")){
	    system("rm -f $chemin_dest".'limits_of_promot_'.$ID_bact."$end_file".".txt");
	}
	open(FSORLIMITS_PROMOT,">> $chemin_dest".'limits_of_promot_'.$ID_bact."$end_file".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".'limits_of_promot_'.$ID_bact."$end_file".".txt\n";
    }
    # open(FSORPROMOT2,">> $chemin_dest".$entete_fasta.'promot_ss_seq'.$ID_bact."$end_file"."_$date".".txt")or die "Impossible to create $chemin_dest".$entete_fasta."promot_ss_seq".$ID_bact."$end_file"."_$date".".txt\n";
    my $long = 0;
   
    my $RBS_upstr = '';
    my $temp_gene_end_prec_for_sens0_limit = 0;
    my $string_limits_sens1 = '';
    my $string_limits_sens0 = '';
    my $prem_sens = undef; # $sens[ $H_ind_T{'gene'} ][0];
    # pour chaque limite inferieure de gene...
    for my $cpt_deb_CDS(0..$#{$beg[ $H_ind_T{'gene'} ]})
    {

	$prem_sens = $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS];

	# ...nous recuperons tous les RBS associes
# 	if(defined $imp_linked_to[$cpt_deb_CDS]){
# 	    # Pour cela, nous listons tous les 'imp' associes a ce gene
# 	    for my $nri(0..$#{ $imp_linked_to[$cpt_deb_CDS] }){
# 		# sens du (des) RBS A VERIFIER !!!
		
# 		for my $nb_sens(0..$#{ $sens[ $H_ind_T{'imp'} ][$nri] }){
# 		    # nous ne devons garder QUE ceux etant des RBS ET DANS LE MEME SENS
# 		    if($sens[$H_ind_T{'gene'}][$cpt_deb_CDS] == $sens[$H_ind_T{'imp'}][$nri][$nb_sens] )
# 		    {
# 			($imp_type[ $imp_linked_to[$cpt_deb_CDS][$nri][1] ] eq 'RBS') and $RBS_upstr .= " beg $beg[ $H_ind_T{'imp'} ][$nri][  $imp_linked_to[$cpt_deb_CDS][0] ], end $end[ $H_ind_T{'imp'} ][$nri][  $imp_linked_to[$cpt_deb_CDS][0] ]";
# 		    }
# 	        }
# 	    }
# 	}


	($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
	my ($gene_beg, $gene_end) = &get_limits($cpt_deb_CDS);

	my ($gene_beg_prec, $gene_end_prec);



# ?????

	    my ($gene_beg_suiv, $gene_end_suiv);

	    if($prem_sens == -1)
	    {
		($not_verbose)or print "\n>Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n"; # , upstr RBS:$RBS_upstr
		
		# si l'on n'a pas affaire au dernier CDS # ??
		if($cpt_deb_CDS != $#{$beg[ $H_ind_T{'gene'} ]})
		{	
		    ($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
		    ($gene_beg_suiv, $gene_end_suiv) = &get_limits($cpt_deb_CDS+1);
		    # longueur = difference entre fin du CDS courant et debut du CDS precedent + la marge (quelques nucleotides apres) 
		    $long = $gene_beg_suiv - $gene_end-1 -$dist_fin_up_seq_deb_trad;
		    print "$prog_tag longueur vaut ".$gene_beg_suiv."-".$gene_end."-1 - $dist_fin_up_seq_deb_trad\n";
		    
		}
		else
		{
		    # longueur = longueur de la sequence globale - limite de fin du gene - la marge souhaitee
		    $long = $length_sequence -$end[ $H_ind_T{'gene'} ][0]-1 -$dist_fin_up_seq_deb_trad;
		    
		}
		
		# sequence inexistante... on ne met que l'en-tete (si l'utilsateur le souhaite)
		if($long < 1 -$dist_fin_up_seq_deb_trad)
		{
		    print FSORPROMOT ">Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n\n"; # , upstr RBS:$RBS_upstr
    
		  #  print FSORPROMOT2 "No intergenic region corresponding to\n";
		    ($not_verbose)or print "$prog_tag No intergenic region corresponding to $gene[$cpt_deb_CDS]\n";
		    ($edit_limits_of_promoters)and $string_limits_sens1 .= "-1,-1;";
		}
		# sequence trop courte par rapport a la limite basse definie (longueurmini)
		elsif($long < $longueurmini)
		{
		    print FSORPROMOT ">Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n\n"; #  upstr RBS:$RBS_upstr    
		    # print FSORPROMOT2 "Intergenic region to short (mini length:$longueurmini, distance from translation start:$dist_fin_up_seq_deb_trad)\n";
		    ($not_verbose)or print "Intergenic region to short (mini length:$longueurmini, distance from translation start:$dist_fin_up_seq_deb_trad)\n";
		    print "$prog_tag longueur $long inferieure a longueur_mini $longueurmini\n";
		    ($edit_limits_of_promoters) and $string_limits_sens1 .= "-1,-1;";
		}
		else
		{
		    if(!defined $gene[$cpt_deb_CDS]){ print "$prog_tag gene not defined for $cpt_deb_CDS\n";}

		    ($not_verbose)or print "$prog_tag on imprime la partie de la seq commencant a la position $gene_end+$dist_fin_up_seq_deb_trad =".($gene_end+$dist_fin_up_seq_deb_trad)." et de longueur $long\n";
		    ($not_verbose)or print &print50($fasta,&complement_seq(substr($sequence,$gene_end+$dist_fin_up_seq_deb_trad, $long)));		    		    
		    
		    print FSORPROMOT ">Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n"; # , upstr RBS:$RBS_upstr    
		    print FSORPROMOT &print50($fasta,&complement_seq(substr($sequence,$gene_end+$dist_fin_up_seq_deb_trad, $long)));
		    ($edit_limits_of_promoters) and $string_limits_sens1 .= $gene_end.",".($gene_end + $long-1).";"
		}
	    }
	    # traitement des sequences en sens 1
	    elsif($prem_sens == 1)
	    {
		($not_verbose)or print "\n>Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n"; # , upstr RBS:$RBS_upstr
		
		my $debut = undef;
		
		# cas general, nous ne traitons pas le premier gene
		if($cpt_deb_CDS)
		{
		    # ??
		    ($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
		    ($gene_beg_prec, $gene_end_prec) = &get_limits($cpt_deb_CDS-1);
		    if($gene_end_prec > $temp_gene_end_prec_for_sens0_limit)
		    {
			$temp_gene_end_prec_for_sens0_limit = $gene_end_prec;		     
		    }
		    else
		    {
			print "$prog_tag LIMITE: pour le gene $gene[$cpt_deb_CDS] la fin du gene précedent $gene_end_prec , est inférieur à la fin du gene d'avant, $temp_gene_end_prec_for_sens0_limit\n";
		    }
		    $debut = $temp_gene_end_prec_for_sens0_limit;
		    $long = $gene_beg - $temp_gene_end_prec_for_sens0_limit - 1 -$dist_fin_up_seq_deb_trad;
		    print "$prog_tag on memorise pour debut $debut et pour longueur $long ($gene_beg - $temp_gene_end_prec_for_sens0_limit - 1 -$dist_fin_up_seq_deb_trad) pour $gene[$cpt_deb_CDS]\n";
		    
		}
		# traitement du premier gene du genome
		else
		{
		    $debut = 0;
		    $long = $beg[ $H_ind_T{'gene'} ][0] - 1-$dist_fin_up_seq_deb_trad;
		    print "$prog_tag on memorise pour debut $debut et pour longueur $long ($beg[ $H_ind_T{'gene'} ][0] - 1 - $dist_fin_up_seq_deb_trad)\n";
		    # print "long $long vaut $beg[ $H_ind_T{'gene'} ][0] - 1 - $dist_fin_up_seq_deb_trad\n";
		    
		}
		
		# sequence inexistante
		if($long < 1-$dist_fin_up_seq_deb_trad )
		{
		    print FSORPROMOT ">Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n\n";  # , upstr RBS:$RBS_upstr
		    # print FSORPROMOT2 "No intergenic region corresponding to\n";
		    ($not_verbose)or print "$prog_tag No intergenic region corresponding to\n";
		    ($edit_limits_of_promoters) and $string_limits_sens0 .= "-1,-1;";
		}
		# sequence trop courte
		elsif($long < $longueurmini)
		{
		    print FSORPROMOT ">Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n\n"; # , upstr RBS:$RBS_upstr    
		    # print FSORPROMOT2 "Intergenic region to short (mini length:$longueurmini, distance from translation start:$dist_fin_up_seq_deb_trad)\n";
		    ($not_verbose)or print "$prog_tag Intergenic region to short ($long)(mini length:$longueurmini, distance from translation start:$dist_fin_up_seq_deb_trad)\n";
		    ($edit_limits_of_promoters) and $string_limits_sens0 .= "-1,-1;";
		}
		# cas general
		else
		{
		    print FSORPROMOT ">Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n"; #  upstr RBS:$RBS_upstr  
		    print "$prog_tag Debut $debut long $long ".length($sequence)."\n";
		    if($debut + $long >= length($sequence)){
			print "$prog_tag pb gene $gene[$cpt_deb_CDS] \n";
			exit;
		    }
		    print FSORPROMOT &print50($fasta,substr($sequence,$debut, $long));	    
		    ($not_verbose)or print &print50($fasta,substr($sequence,$debut, $long));
		    print "$prog_tag on memorise pour debut $debut et pour longueur $long\n";
		    # die "on extrait de $debut a ".($beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] -1)." pour $gene[$cpt_deb_CDS] sens $prem_sens\n";
		    ($edit_limits_of_promoters) and $string_limits_sens0 .= $debut.",".($debut + $long-1).";";
		}
	    }
	    # sens non connu
	    else{
		print FSORPROMOT ">Promot_for:"."$gene[$cpt_deb_CDS] sens $prem_sens (unknown) deb_trad ".(-$dist_fin_up_seq_deb_trad)."\n\n"; #  upstr RBS:$RBS_upstr
		print "$prog_tag undefined sens for $prem_sens line ".__LINE__."\n";
	    }

    }
    
    close FSORPROMOT;
    # close FSORPROMOT2;
    print "$prog_tag $chemin_dest".$entete_fasta.'promot_'.$ID_bact."$end_file".".txt file created\n";
    # print "$chemin_dest".$entete_fasta.'promot_ss_seq'.$ID_bact."$end_file"."_$date".".txt file created\n";
    if($edit_limits_of_promoters)
    {
	print FSORLIMITS_PROMOT $string_limits_sens0."\n".$string_limits_sens1."\n";
	close FSORLIMITS_PROMOT;
    }
    
}


if($edit_file_for_all_interCDS)
{
#     if(($bool_premier_passage) and (-e "$chemin_dest".$entete_fasta.'seq_interCDS_'.$ID_bact."$end_file"."_$date".".txt") )
#     {
# 	die "File $chemin_dest".$entete_fasta.'seq_interCDS_'.$ID_bact."$end_file"."_$date".".txt already exists, change the destination folder \n";
#     }

    if($bool_premier_passage and (-e "$chemin_dest".$entete_fasta.'seq_interCDS_'.$ID_bact."$end_file"."_$date".".txt")){ 
	system("rm -f $chemin_dest".$entete_fasta.'seq_interCDS_'.$ID_bact."$end_file"."_$date".".txt");
    }
    # pour toutes les sequences interCDS
    open(FSORENSINTERCDS,">> $chemin_dest".$entete_fasta.'seq_interCDS_'.$ID_bact."$end_file"."_$date".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".$entete_fasta."seq_interCDS_".$ID_bact."$end_file"."_$date".".txt\n";
    
    if($beg[ $H_ind_T{'gene'} ][0] != 1)
    {
	$cpt_seq++;
	print FSORENSINTERCDS ">InterCDS_sequence_$cpt_seq\n";
	print FSORENSINTERCDS &print50($fasta,substr($sequence,0,$beg[ $H_ind_T{'gene'} ][0]-1));
	($not_verbose)or print "\n>InterCDS sequence $cpt_seq 1\n";
	($not_verbose)or print &print50($fasta,substr($sequence,0,$beg[ $H_ind_T{'gene'} ][0]-1));
    }

    for my $chaqueCDS(0..($#{$beg[ $H_ind_T{'gene'} ]}-1))
    {
	($not_verbose)or print "\ndeb_CDS fin_CDS\n$beg[ $H_ind_T{'gene'} ][$chaqueCDS+1] $end[ $H_ind_T{'gene'} ][$chaqueCDS]";
	if( (my $long = $beg[ $H_ind_T{'gene'} ][$chaqueCDS+1]-$end[ $H_ind_T{'gene'} ][$chaqueCDS]-1)>0 )
	{
	    $cpt_seq++;
	    print FSORENSINTERCDS ">InterCDS_sequence_$cpt_seq\n";
	    print FSORENSINTERCDS &print50($fasta,substr($sequence,$end[ $H_ind_T{'gene'} ][$chaqueCDS], $long));
	    ($not_verbose)or print "\n>InterCDS sequence $cpt_seq\n";
	    ($not_verbose)or print &print50($fasta,substr($sequence,$end[ $H_ind_T{'gene'} ][$chaqueCDS], $long));
	}
	else
	{
	    ($not_verbose)or print "\n$prog_tag Problem, upstream limit of interCDS  higher than downstream limit!\nCSD must overlap!\n";
	    # print FPB "Problem, upstream limit of interCDS  higher than downstream limit!\nCSD must overlap!\n";
	}
    }
    
    if($length_sequence-$end[ $H_ind_T{'gene'} ][$#{$end[ $H_ind_T{'gene'} ]}] != 0)
    {
	$cpt_seq++;
	print FSORENSINTERCDS ">InterCDS_sequence_$cpt_seq\n";
	print FSORENSINTERCDS &print50($fasta,substr($sequence,$end[ $H_ind_T{'gene'} ][$#{$end[ $H_ind_T{'gene'} ]}],$length_sequence - $end[ $H_ind_T{'gene'} ][$#{$end[ $H_ind_T{'gene'} ]}]));
	
	($not_verbose)or print "\n>InterCDS sequence $cpt_seq\n";
	($not_verbose)or print &print50($fasta,substr($sequence,$end[ $H_ind_T{'gene'} ][$#{$end[ $H_ind_T{'gene'} ]}],$length_sequence - $end[ $H_ind_T{'gene'} ][$#{$end[ $H_ind_T{'gene'} ]}]));
    }
    close FSORENSINTERCDS;
    print "$prog_tag $chemin_dest".$entete_fasta.'seq_interCDS_'.$ID_bact."$end_file"."_$date".".txt file created\n";
    $cpt_seq = 0;
}

# if($edit_file_for_all_interCDS_between_conv)
# {
    # pour les sequences interCDS entre genes convergents
   # open(FSORCONV,"> $chemin_dest".$entete_fasta.'seq_interCDS_entre_conv_'.$ID_bact."$end_file"."_$date".".txt")or die "Impossible to create $chemin_dest".$entete_fasta."seq_interCDS_entre_conv_".$ID_bact."$end_file"."_$date".".txt\n";
    
    
  #  for my $chaqueCDS(0..$#deb_interCDS_entre_conv)
  #  {
	# my $long = $fin_interCDS_entre_conv[$chaqueCDS]-$deb_interCDS_entre_conv[$chaqueCDS];
	# $cpt_seq++;
	# print FSORCONV ">InterCDS_between_convergent_genes_$cpt_seq\n";
	# print FSORCONV &print50($fasta,substr($sequence,$deb_interCDS_entre_conv[$chaqueCDS], $long));
	# ($not_verbose)or print "\n>InterCDS between convergent genes $cpt_seq\n";
	# ($not_verbose)or print &print50($fasta,substr($sequence,$deb_interCDS_entre_conv[$chaqueCDS], $long));
    # }
    # close FSORCONV;
    # print "$chemin_dest".'seq_interCDS_entre_conv_'.$ID_bact."$end_file"."_$date".".txt file created\n";
    # $cpt_seq = 0;
# }

# if($edit_file_for_all_interCDS_between_div)
# {
    # pour les sequences interCDS entre genes divergents
    # open(FSORDIV,"> $chemin_dest".$entete_fasta.'seq_interCDS_entre_div_'.$ID_bact."$end_file"."_$date".".txt")or die "Impossible to create $chemin_dest".$entete_fasta."seq_interCDS_entre_div_".$ID_bact."$end_file"."_$date".".txt\n";
    
    # for my $chaqueCDS(0..$#deb_interCDS_entre_div)
    # {
	# my $long = $fin_interCDS_entre_div[$chaqueCDS]-$deb_interCDS_entre_div[$chaqueCDS];
	# $cpt_seq++;
	# print FSORDIV ">InterCDS_between_divergent_genes_$cpt_seq\n";
	# print FSORDIV &print50($fasta,substr($sequence,$deb_interCDS_entre_div[$chaqueCDS], $long));
	# ($not_verbose)or print "\n>InterCDS between divergent genes $cpt_seq\n";
#	($not_verbose)or print &print50($fasta,substr($sequence,$deb_interCDS_entre_div[$chaqueCDS], $long));
  #  }
  #  close FSORDIV;
  #  print "$chemin_dest".'seq_interCDS_entre_div_'.$ID_bact."$end_file"."_$date".".txt file created\n";
  #  $cpt_seq = 0;
# }

#if($edit_file_for_all_CDS_interCDS_between_conv)
#{
    #pour les sequence CDS+interCDS entre genes conv
 #   open(FSORCDS_INTERCDSCONV,"> $chemin_dest".$entete_fasta.'seq_CDS_et_interCDS_entre_conv_'.$ID_bact."$end_file"."_$date".".txt")or die "Impossible to create $chemin_dest".$entete_fasta."seq_CDS_et_interCDS_entre_conv_".$ID_bact."$end_file"."_$date".".txt\n";
      
#    for my $chaqueCDS(0..$#deb_CDS_interCDS_entre_conv)
#    {
#        my $long = $fin_CDS_interCDS_entre_conv[$chaqueCDS]-$deb_CDS_interCDS_entre_conv[$chaqueCDS]+1;
#        $cpt_seq++;
#        print FSORCDS_INTERCDSCONV ">CDS&interCDS_between_convergent_genes_$cpt_seq\n";
#	 print FSORCDS_INTERCDSCONV &print50($fasta,substr($sequence,$deb_CDS_interCDS_entre_conv[$chaqueCDS]-1, $long));
#	 ($not_verbose)or print "\n>CDS&interCDS between convergent genes $cpt_seq\n";
#	 ($not_verbose)or print &print50($fasta,substr($sequence,$deb_CDS_interCDS_entre_conv[$chaqueCDS]-1, $long));	
#    }
#    close FSORCDS_INTERCDSCONV;
#    print "$chemin_dest".$entete_fasta.'seq_CDS_et_interCDS_entre_conv_'.$ID_bact."$end_file"."_$date".".txt file created\n";
#    $cpt_seq = 0;
#}

# OK
if($edit_file_for_all_annotations)
{
    my $bool_with_limits = 1;
    # pour toutes les annotations

#     if(($bool_premier_passage) and (-e "$chemin_dest".'annot_'.$ID_bact.$end_file.".txt") )
#     {
# 	die "File $chemin_dest".'annot_'.$ID_bact.$end_file.".txt already exists, change the destination folder \n";
#     }
#     

    if($bool_premier_passage and (-e "$chemin_dest".'annot_'.$ID_bact.$end_file.".txt")){ 
	system("rm -f $chemin_dest".'annot_'.$ID_bact.$end_file.".txt");
    }

    # pour toutes les sequences CDS (donc leurs annotations)
    open(FSOREACHANNOT,">> $chemin_dest".'annot_'.$ID_bact.$end_file.".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest"."annot_".$ID_bact.$end_file.".txt\n";
    for my $numg(0..$#gene)
    {
	# my $long = $end[ $H_ind_T{'gene'} ][$numg] - $beg[ $H_ind_T{'gene'} ][$numg] +1;
	print FSOREACHANNOT ">annot_for_$gene[$numg]";    
	# (defined $prot_name[$numg])   and print FSOREACHANNOT ", prot_name $prot_name[$numg], Nr_acc_prot $Nr_acc_prot[$numg], ";	
	print FSOREACHANNOT ", prot_name $prot_name[$numg]\n";
	# print FSOREACHANNOT "\n";    
	# ($not_verbose)or print "\n>annot for $gene[$numg]";
	
# 	if(defined $imp_linked_to[$numg]){
# 	    # print "DEFINI\n";1
# 	    print FSOREACHANNOT "RELATED ANNOTATIONS:\n";
# 	    for my $nri(0..$#{ $imp_linked_to[$numg] }){
# 		# print "nri $nri\n";
# 		# print FSOREACHANNOT "imp: $imp_type_kind[ $imp_type[ $imp_linked_to[$nri][1] ]] , sens $sens[ $H_ind_T{'imp'} ][ $imp_linked_to[$nri][1] ][ $imp_linked_to[$nri][0] ], beg $beg[ $H_ind_T{'imp'} ][ $imp_linked_to[$nri][1] ][ $imp_linked_to[$nri][0] ], end $end[ $H_ind_T{'imp'} ][ $imp_linked_to[$nri][1] ][ $imp_linked_to[$nri][0] ], ";
# 		# print "imp_linked_to nri $nri 1 $imp_linked_to[$numg][$nri][1]\n";    
# 		# print "imp_linked_to nri $nri 0 $imp_linked_to[$numg][$nri][0]\n";    
# 		# print "imp: imp_type $imp_type[ $imp_linked_to[$numg][$nri][1] ]\n";
# 		# print "sens sens $sens[ $H_ind_T{'imp'} ][ $imp_linked_to[$numg][$nri][1] ][ $imp_linked_to[$numg][$nri][0] ]\n";
# 		print FSOREACHANNOT "$imp_type_kind[ $imp_type[ $imp_linked_to[$numg][$nri][1] ]] , imp_sens $sens[ $H_ind_T{'imp'} ][ $imp_linked_to[$numg][$nri][1] ][ $imp_linked_to[$numg][$nri][0] ], $beg[ $H_ind_T{'imp'} ][ $imp_linked_to[$numg][$nri][1] ][ $imp_linked_to[$numg][$nri][0] ]..$end[ $H_ind_T{'imp'} ][ $imp_linked_to[$numg][$nri][1] ][ $imp_linked_to[$numg][$nri][0] ]";
# 		defined($comment[ $H_ind_T{'imp'} ][$nri] ) and do 
# 		{ 
# #		defined($comment[ $H_ind_T{'imp'} ][ $imp_linked_to[$numg][$nri] ]) and do { 
# 		    foreach(@{ $comment[ $H_ind_T{'imp'} ][$nri] }){
# 			print FSOREACHANNOT ", comment_imp: $_"; 
# 		    }
# 		};	
	
# 		print FSOREACHANNOT "\n";
# 		# }
# 	    }
# 	    # exit;
# 	}
# 	else{ print "imp_linked_to of $gene[$numg]  NOT DEF\n"; }
# ?????

	    (defined $sens[ $H_ind_T{'gene'} ][$numg])and do { ($sens[ $H_ind_T{'gene'} ][$numg] == 3)and print "$prog_tag PB: UNKNOWN SENS for $gene[$numg]\n"; 
		  print FSOREACHANNOT "gene_sens $sens[ $H_ind_T{'gene'} ][$numg], "; 
								    };
	    ($bool_with_limits)and print FSOREACHANNOT "$beg[ $H_ind_T{'gene'} ][$numg]..$end[ $H_ind_T{'gene'} ][$numg], ";

	(defined $Nr_acc_prot[$numg])   and print FSOREACHANNOT "Nr_acc_prot $Nr_acc_prot[$numg], ";
#	(defined $seqdesc_title[ $H_ind_T{'prot'} ][$numg])and print FSOREACHANNOT "seqdesc_title $seqdesc_title[ $H_ind_T{'prot'} ][$numg], ";
	

		print FSOREACHANNOT "$beg[ $H_ind_T{'prot'} ][$numg]..$end[ $H_ind_T{'prot'} ][$numg], ";
		(defined $sens[ $H_ind_T{'prot'} ][$numg])and print FSOREACHANNOT "prot_sens $sens[ $H_ind_T{'prot'} ][$numg], ";


	# print 'prot '.substr($prot[$numg],0,7)."\n"; 
	if (defined($comment[ $H_ind_T{'gene'} ][$numg])){ 
		print FSOREACHANNOT "\ngene comment :$comment[ $H_ind_T{'gene'} ][$numg]; ";
	}
	
	if (defined($comment[ $H_ind_T{'prot'} ][$numg])){
		print FSOREACHANNOT "\nprot comment :$comment[ $H_ind_T{'prot'} ][$numg]; ";
	}
	
# 	if((defined $comment_cdregion[$numg]) and ($comment_cdregion[$numg] ne ''))
# 	{
# 	    print FSOREACHANNOT "\ncomment_cdregion : $comment_cdregion[$numg]";
# 	    print " comment cdregion $comment_cdregion[$numg] pour $numg $gene[$numg]\n";
# 	}

# 	if(defined $region_linked_to[$numg])
# 	{
	   
# 	    print FSOREACHANNOT "\nRegion(s) :";	    
# 	    for my $nri(0..$#{ $region_linked_to[$numg] })	   
# 	    {	
# 		print FSOREACHANNOT "$region[$region_linked_to[$numg][$nri]],";					
# 	    }
	   
# 	}

	print FSOREACHANNOT "\n";
	# ($not_verbose)or print ;
    }
    close FSOREACHANNOT;
    print $prog_tag.' '.$chemin_dest.'annot_'.$ID_bact.$end_file.".txt file created\n";
#    $cpt_seq = 0;
}


if($edit_file_for_upstream_seq)
{
# !!!!!!!ATTENTION: cas deb trop courts par rapport a longueur demandee non encore teste!!!!!!!!!!!!!!!!!! (decallage possible d'une lettre?)
#     if(($bool_premier_passage) and (-e "$chemin_dest".$entete_fasta."upstr_seq_".(-$dist_deb_up_seq_deb_trad)."_".(-$dist_fin_up_seq_deb_trad).'_'.$ID_bact."$end_file"."_$date".".txt") )
#     {
# 	die "File $chemin_dest".$entete_fasta."upstr_seq_".(-$dist_deb_up_seq_deb_trad)."_".(-$dist_fin_up_seq_deb_trad).'_'.$ID_bact."$end_file"."_$date".".txt already exists, change the destination folder \n";
#     }
#     
    if($bool_premier_passage and (-e "$chemin_dest".$entete_fasta."upstr_seq_".(-$dist_deb_up_seq_deb_trad)."_".(-$dist_fin_up_seq_deb_trad).'_'.$ID_bact."$end_file".".txt")){ 
	system("rm -f $chemin_dest".$entete_fasta."upstr_seq_".(-$dist_deb_up_seq_deb_trad)."_".(-$dist_fin_up_seq_deb_trad).'_'.$ID_bact."$end_file".".txt");
    }

    # pour toutes les sequences amont (chevauchements possibles)
    open(FSORUPSEQ,">> $chemin_dest".$entete_fasta."upstr_seq_".(-$dist_deb_up_seq_deb_trad)."_".(-$dist_fin_up_seq_deb_trad).'_'.$ID_bact."$end_file".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".$entete_fasta."upstr_seq_".(-$dist_deb_up_seq_deb_trad)."_".(-$dist_fin_up_seq_deb_trad).'_'.$ID_bact."$end_file".".txt\n";
   
    for my $cpt_deb_CDS(0..$#{$beg[ $H_ind_T{'gene'} ]})
    {
	my $RBS_upstr = 0;

# 	if(defined $imp_linked_to[$cpt_deb_CDS]){
# 	    for my $nri(0..$#{ $imp_linked_to[$cpt_deb_CDS] }){
# 		# for my $nb_sens(0..$#{ $sens[ $H_ind_T{'imp'} ][$nri] }){
# 		($imp_type_kind[$imp_type[ $imp_linked_to[$cpt_deb_CDS][$nri][1] ]] eq 'RBS') and $RBS_upstr .= " beg $beg[ $H_ind_T{'imp'} ][ $imp_linked_to[$cpt_deb_CDS][$nri][1] ][ $imp_linked_to[$cpt_deb_CDS][$nri][0] ], end $end[ $H_ind_T{'imp'} ][ $imp_linked_to[$cpt_deb_CDS][$nri][1] ][ $imp_linked_to[$cpt_deb_CDS][$nri][0] ]";
# 		# }
# 	    }
# 	}
	

	    
	    #if(exists $H_RBS{$cpt_deb_CDS}){ $RBS_upstr = $H_RBS{$cpt_deb_CDS};}
	    # traitement des sequences en sens -1
	    
	    if($sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == -1 or $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == 2)
	    {	
		($not_verbose)or print "\n>Promot_for:"."$gene[$cpt_deb_CDS] sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] deb_trad ".(-$dist_fin_up_seq_deb_trad).", upstr RBS:$RBS_upstr\n";  
		#si ce n'est pas le dernier gene
		if($length_sequence - $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS] >= $dist_deb_up_seq_deb_trad)
		{
		    print FSORUPSEQ ">Promot_for:"."$gene[$cpt_deb_CDS] sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] deb_trad ".(-$dist_fin_up_seq_deb_trad).", upstr RBS:$RBS_upstr\n";  
		    print FSORUPSEQ &print50($fasta,&complement_seq(substr($sequence,$end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]+$dist_fin_up_seq_deb_trad, $length_upstream)));
		    ($not_verbose)or print &print50($fasta,&complement_seq(substr($sequence,$end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]+$dist_fin_up_seq_deb_trad, $length_upstream)));   
		}
		elsif((my $length_end = $length_sequence - $dist_fin_up_seq_deb_trad - $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]) >= $longueurmini)
		{
		    print FSORUPSEQ ">Promot_for:"."$gene[$cpt_deb_CDS] sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS]deb_trad ".(-$dist_fin_up_seq_deb_trad).", upstr RBS:$RBS_upstr\n";  
		    print FSORUPSEQ &print50($fasta,&complement_seq(substr($sequence,$end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]+$dist_fin_up_seq_deb_trad, $length_end)));
		    ($not_verbose)or print &print50($fasta,&complement_seq(substr($sequence,$end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]+$dist_fin_up_seq_deb_trad, $length_end)));   
		}
		else
		{
		    print "$prog_tag Upstream sequence of $gene[$cpt_deb_CDS] gene to short to be used\n";
		}
	    }
	    elsif($sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == 1 or $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == 2)
	    {	 
		# traitement des sequences en sens 1  
		
		($not_verbose)or print "\n>Promot_for:"."$gene[$cpt_deb_CDS] sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] upstr RBS:$RBS_upstr\n";  
		if( $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] > $dist_deb_up_seq_deb_trad)
		{
		    print FSORUPSEQ ">Promot_for:"."$gene[$cpt_deb_CDS] sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] deb_trad ".(-$dist_fin_up_seq_deb_trad).", upstr RBS:$RBS_upstr\n";  
		    print FSORUPSEQ &print50($fasta,substr($sequence,$beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - $dist_deb_up_seq_deb_trad -1, $length_upstream));
		    ($not_verbose)or print &print50($fasta,substr($sequence,$beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - $dist_deb_up_seq_deb_trad -1, $length_upstream));
		}
		elsif($beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] > $dist_fin_up_seq_deb_trad + $longueurmini)
		{
		    print FSORUPSEQ ">Promot_for:"."$gene[$cpt_deb_CDS] sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] deb_trad ".(-$dist_fin_up_seq_deb_trad).", upstr RBS:$RBS_upstr\n";  
		    print FSORUPSEQ &print50($fasta,substr($sequence, 0, $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - $dist_fin_up_seq_deb_trad));
		    ($not_verbose)or print &print50($fasta,substr($sequence, 0, $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - $dist_fin_up_seq_deb_trad));
		    print "$prog_tag Verifier la longueur de la sequence amont du premier gene\n";
		}
		else
		{
		    print "$prog_tag Upstream sequence of $gene[$cpt_deb_CDS] gene to short to be used\n";
		}
	    }
	    else{
		die "$prog_tag [Error] case of unknown sens $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] we have to treat line ".__LINE__."\n";
	    }

    }
    close FSORUPSEQ;

    print "$prog_tag $chemin_dest".$entete_fasta."upstr_seq_".(-$dist_deb_up_seq_deb_trad)."_".(-$dist_fin_up_seq_deb_trad).'_'.$ID_bact."$end_file".".txt file created\n";  
#    $cpt_seq = 0;
}

# edition des sequences proteiques
if($edit_file_for_protein_seq)
{
#     if(($bool_premier_passage) and (-e "$chemin_dest".$entete_fasta.'prot_seq_'.$ID_bact."$end_file".".txt") )
#     {
# 	die "File $chemin_dest".$entete_fasta.'prot_seq_'.$ID_bact."$end_file".".txt already exists, change the destination folder \n";
#     }
    if($bool_premier_passage and (-e "$chemin_dest".$entete_fasta.'prot_seq_'.$ID_bact."$end_file".".txt")){ 
	system("rm -f $chemin_dest".$entete_fasta.'prot_seq_'.$ID_bact."$end_file".".txt");
    }
    open(FSORPROTSEQ,">> $chemin_dest".$entete_fasta.'prot_seq_'.$ID_bact."$end_file".".txt")or die "$prog_tag [Error] Impossible to create $chemin_dest".$entete_fasta."prot_seq_".$ID_bact."$end_file".".txt\n";
    # print "idgene $#gene seqprot $#prot\n";
    for my $prot(0..$#gene)
    { 
	if(! defined $prot[$prot])
	{
	    print "$prog_tag seqprot de prot non defini pour prot egal $prot correspondant a $gene[$prot] $prot $#prot\n";
	}
	if(! defined $gene[$prot])
	{
	    print "$prog_tag gene de geneprot de prot non defini pour prot egal $prot\n";
	}
	print FSORPROTSEQ ">$gene[$prot] ".$prot_name[$prot]." Nr_acc_prot $Nr_acc_prot[$prot]\n".&print50($fasta,$prot[$prot]);
	if(not defined $prot[$prot]){ die "$prog_tag [Error] prot not defined for gene $gene[$prot]\n"; }
    } 
    close FSORPROTSEQ;
    print "$prog_tag $chemin_dest".$entete_fasta.'prot_seq_'.$ID_bact."$end_file"."_$date".".txt file created\n";  
}


if($edit_file_for_300_0_seq or $edit_limits_of_overlapped_seq_sens01 or $edit_annotations_of_limits_overlapped){

    # edition des sequences -300 0 en sens 1 puis -1 chevauchantes (pour stat)
    my ($ind0,$ind1) = (undef, undef);            # variables qui vont stocker les indices des premiers gènes respectivement en sens 1 et en sens -1 (pour l'initialisation des variables lors de l'extraction des séquences -300 0 chevauchantes (blocs pour stats)
    
    for my $g(0..$#{ $beg[ $H_ind_T{'gene'} ] }){

	    if((defined $ind0)and(defined $ind1)){ print "$prog_tag ind0 $ind0, ind1 $ind1 dans boucle, line ".__LINE__."\n";last; } 
	    if   (($sens[ $H_ind_T{'gene'} ][$g] == -1)and(not defined $ind1)){ $ind1 = $g; print "$prog_tag ind1 $ind1, line ".__LINE__."\n"; next; }
	    if(($sens[ $H_ind_T{'gene'} ][$g] == 1)and(not defined $ind0)){ $ind0 = $g; print "$prog_tag ind0 $ind0, line ".__LINE__."\n"; next; }
# 	    print "sens exotique $sens[ $H_ind_T{'gene'} ][$g]\n"; die "error line ".__LINE__."\n"; 

    }
    #print "ind1 $ind1 ind0 $ind0 \n";
 

    # @current_pos_in_files = ();

    if($bool_premier_passage){
      
      my @diff_file_title = (
			     '_overlapped_seq_sens0_',
			     '_overlapped_seq_sens1_',
			     '_limits_of_overlapped_seq_sens01_'
			    );
      my @diff_merged_file_title = (
				    '_overlapped_seq_sens01_',
				   );      

      if($edit_annotations_of_limits_overlapped){
	push @diff_file_title, (
				'_annot_overlapped_seq_sens0_',
				'_annot_overlapped_seq_sens1_',
			       );
	push @diff_merged_file_title, (
				      '_annot_overlapped_seq_sens01_'
				     );
      }
      
      $last_ind_1sens_file = $#diff_file_title;

      foreach(@diff_file_title, @diff_merged_file_title){
	push @tmp_file_name, "$chemin_dest".$entete_fasta.(-$dist_deb_up_seq_deb_trad).'_'.(-$dist_fin_up_seq_deb_trad).$_.$ID_bact."$end_file".".txt";
	print "$prog_tag test for $tmp_file_name[$#tmp_file_name] file name\n";
	if(-e $tmp_file_name[$#tmp_file_name]){
	  print "$prog_tag File $tmp_file_name[$#tmp_file_name] already exists, change the destination folder, we remove\n";
	  system("rm -f $tmp_file_name[$#tmp_file_name]");
	}
      }
    }
    
    open(FSOR300_0SEQ ,">> $tmp_file_name[0]")or die "$prog_tag [Error] Impossible to create $tmp_file_name[0] line ".__LINE__."\n"; # _overlapped_seq_sens0_
    open(FSOR300_1SEQ ,">> $tmp_file_name[1]")or die "$prog_tag [Error] Impossible to create $tmp_file_name[1] line ".__LINE__."\n"; # _overlapped_seq_sens1_
    open(FSORLIMSEQ   ,">> $tmp_file_name[2]")or die "$prog_tag [Error] Impossible to create $tmp_file_name[2] line ".__LINE__."\n"; # _limits_of_overlapped_seq_sens01_

    if($edit_annotations_of_limits_overlapped){

      open(FSORANNOT_OVERLAPPED_SENS0 ,">> $tmp_file_name[3]")or die "$prog_tag [Error] Impossible to create $tmp_file_name[3] line ".__LINE__."\n"; # _annot_overlapped_seq_sens0_
      open(FSORANNOT_OVERLAPPED_SENS1 ,">> $tmp_file_name[4]")or die "$prog_tag [Error] Impossible to create $tmp_file_name[4] line ".__LINE__."\n"; # _annot_overlapped_seq_sens1_
    }

    if(not $not_verbose){ 
       foreach(@tmp_file_name[0..$last_ind_1sens_file]){ print "$prog_tag $_ file opened\n"; } 
    }

    if(! defined $end[ $H_ind_T{'gene'} ][0]){ die "$prog_tag [Error] Problem, no CDS limits found by the programme\n"; }
    # variable used to detect overlaps in the same sens (-1)
    ($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
    
    my ($gene_beg_ind1, $gene_end_ind1) = (undef, undef);
    my ($prev_beg1,$prev_end1,$ID_1) = (undef, undef);
    if(defined $ind1){
      
      ($gene_beg_ind1, $gene_end_ind1) = &get_limits($ind1); 
      $prev_beg1 = $gene_end_ind1 + $dist_fin_up_seq_deb_trad; 
      $prev_end1 = $gene_end_ind1 + $dist_deb_up_seq_deb_trad -1; # prev_end1 record the end of the previous upstream seq (the begin in biological terminology)
      $ID_1 = '';
    }
    
    # variable used to detect overlaps in the same sens (1)
    ($not_verbose) or print "$prog_tag Appel get_limits line ".__LINE__."\n";
    my ($gene_beg_ind0, $gene_end_ind0) = (undef, undef);
    my ($prev_beg0,$prev_end0,$ID_0) = (undef, undef);

    if(defined $ind0){

      ($gene_beg_ind0, $gene_end_ind0) = &get_limits($ind0);
      $prev_end0 = $gene_beg_ind0 - $dist_fin_up_seq_deb_trad - 2;
      $prev_beg0 = $gene_beg_ind0 - $dist_deb_up_seq_deb_trad - 1;
      ($prev_beg0 >= 0)or $prev_beg0 = 0;
      $ID_0 = '';
    }
    
    my $ind_last_letter = $length_sequence - 1;
    my $chaine_sens1 = "";
    my $chaine_sens0 = "";
    my $chaine_annot_sens1 = "";
    my $chaine_annot_sens0 = "";
    my @id_numeric_gene_sens1 = ();
    my @id_numeric_gene_sens0 = ();
    my ($dernier_deb_CDS_sens1,$dernier_deb_CDS_sens0)  = (undef, undef);
    
    # on recherche les dernier gene en sens -1 et en sens 1
  PREMIERE:for(my $cpt_deb_CDS=$#{$beg[ $H_ind_T{'gene'} ]}; $cpt_deb_CDS>=0; $cpt_deb_CDS--){
      

	
	if(( $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == -1) and (not defined $dernier_deb_CDS_sens1)){ 
	  $dernier_deb_CDS_sens1 = $cpt_deb_CDS; 
	}
	elsif(( $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == 1) and (not defined $dernier_deb_CDS_sens0)){
	  $dernier_deb_CDS_sens0 = $cpt_deb_CDS; 
	}
	
	if((defined $dernier_deb_CDS_sens0) and (defined $dernier_deb_CDS_sens1)){
	  last PREMIERE; 
	}

    }

    for my $cpt_deb_CDS(0..$#{$beg[ $H_ind_T{'gene'} ]})
    {

	    my $RBS_upstr = 0;

	    # traitement des sequences en sens -1
	    if( $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == -1 or $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == 2  )
	    {
		if (not defined $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]){ die "$prog_tag [Error] limite end $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS] pas definit\n";}
		
		if( $prev_end1 >= $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]-1+$dist_fin_up_seq_deb_trad )
		{
		    # print "prev_end1 passe de $prev_end1 a ";
		    # we change the limits only if it is an higher one
		    if($prev_end1 < $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS] + $dist_deb_up_seq_deb_trad - 1)
		    {
			$prev_end1 = $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS] + $dist_deb_up_seq_deb_trad - 1;
		    }
		    if($prev_beg1 > $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]  + $dist_fin_up_seq_deb_trad)
		    {
			$prev_beg1 = $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS]  + $dist_fin_up_seq_deb_trad;
		    }
		    # print "$prev_end1 pour $ID_1 ligne 1389\n";
		    ($ID_1 eq '')or $ID_1 .= $separateur;
		    $ID_1 .= $gene[$cpt_deb_CDS];
		    push @id_numeric_gene_sens1, $cpt_deb_CDS;
		    if(($cpt_deb_CDS == $#{$beg[ $H_ind_T{'gene'} ]}) or ($cpt_deb_CDS == $dernier_deb_CDS_sens1))
		    {
			$cpt_seq++;
			($prev_end1 < $ind_last_letter)or $prev_end1 = $ind_last_letter;
			print FSOR300_1SEQ ">contigues_up_seq_sens1_nbr_$cpt_seq dist_deb_trad ".(-$dist_fin_up_seq_deb_trad).": $ID_1".":\n".&print50($fasta,&complement_seq(substr($sequence,$prev_beg1,$prev_end1 - $prev_beg1 + 1)));
			 print "$prog_tag on ecrit deb $prev_beg1, longueur ".($prev_end1 - $prev_beg1 + 1)." qui correspond a $prev_end1 - $prev_beg1 + 1 ligne ".__LINE__." gene : $ID_1\n";
			
			if($edit_limits_of_overlapped_seq_sens01){ $chaine_sens1 .= $prev_beg1.",".$prev_end1.";";}
			if($edit_annotations_of_limits_overlapped){  print FSORANNOT_OVERLAPPED_SENS1 &get_string_annot($ID_1,\@id_numeric_gene_sens1);	}
		    }
		}
		else
		{
		    $cpt_seq++;
		    ($prev_end1 < $ind_last_letter)or $prev_end1 = $ind_last_letter;
		    print FSOR300_1SEQ ">contigues_up_seq_sens1_nbr_$cpt_seq dist_deb_trad ".(-$dist_fin_up_seq_deb_trad).": $ID_1".":\n".&print50($fasta,&complement_seq(substr($sequence,$prev_beg1,$prev_end1 - $prev_beg1 + 1)));
		    print "$prog_tag on ecrit deb $prev_beg1, longueur ".($prev_end1 - $prev_beg1 + 1)." qui correspond a $prev_end1 - $prev_beg1 + 1 ligne ".__LINE__." gene : $ID_1\n";

		    if($edit_limits_of_overlapped_seq_sens01){$chaine_sens1 .= $prev_beg1.",".$prev_end1.";"; }

		    if($edit_annotations_of_limits_overlapped){print FSORANNOT_OVERLAPPED_SENS1 &get_string_annot($ID_1,\@id_numeric_gene_sens1); }

		    @id_numeric_gene_sens1 = ();
		    $ID_1 = "";
		    $ID_1 .= $gene[$cpt_deb_CDS];
		    $prev_beg1 = $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS] + $dist_fin_up_seq_deb_trad;
		    $prev_end1 = $end[ $H_ind_T{'gene'} ][$cpt_deb_CDS] + $dist_deb_up_seq_deb_trad - 1;
		    print "$prog_tag prev_beg1 vaut $prev_beg1, prev_end1 vaut $prev_end1 ligne ".__LINE__."\n"; 
		    push @id_numeric_gene_sens1, $cpt_deb_CDS;	

		    if(($cpt_deb_CDS == $#{$beg[ $H_ind_T{'gene'} ]}) or ($cpt_deb_CDS == $dernier_deb_CDS_sens1))
		    {
			$cpt_seq++;
			($prev_end1 < $ind_last_letter)or $prev_end1 = $ind_last_letter;
			print FSOR300_1SEQ ">contigues_up_seq_sens1_nbr_$cpt_seq dist_deb_trad ".(-$dist_fin_up_seq_deb_trad).": $ID_1".":\n".&print50($fasta,&complement_seq(substr($sequence,$prev_beg1,$prev_end1 - $prev_beg1 + 1)));
			 print "$prog_tag on ecrit deb $prev_beg1, longueur ".($prev_end1 - $prev_beg1 + 1)." qui correspond a $prev_end1 - $prev_beg1 + 1 ligne ".__LINE__." gene : $ID_1\n";
			if($edit_limits_of_overlapped_seq_sens01){$chaine_sens1 .= $prev_beg1.",".$prev_end1.";"; }
			if($edit_annotations_of_limits_overlapped){print FSORANNOT_OVERLAPPED_SENS1 &get_string_annot($ID_1,\@id_numeric_gene_sens1); }
		    }
		    
		}
		
	    }

	    if( $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == 1 or $sens[ $H_ind_T{'gene'} ][$cpt_deb_CDS] == 2 ) 
	    {	 
		# ??ICI A CORRIGER
		# traitement des sequences en sens 1  
		print "$prog_tag $prev_end0 >= $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS]-$dist_fin_up_seq_deb_trad-1?\n";
		if( $prev_end0 >= $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS]-$dist_deb_up_seq_deb_trad-1 )
		{
		    # ??ICI
		    # we change the limit only if it is an higher one
		    if($prev_end0 < $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - 2 - $dist_fin_up_seq_deb_trad)
		    {
			$prev_end0 = $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - 2 - $dist_fin_up_seq_deb_trad;
			print "$prog_tag on met prev_end0 a $prev_end0 line ".__LINE__."\n";
		    }
		    #Dans le sens 1 pas possible
		    if($prev_beg0 > $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - $dist_deb_up_seq_deb_trad - 1)
		    {
			$prev_beg0 = $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - $dist_deb_up_seq_deb_trad - 1;
			print "$prog_tag on met prev_beg0 a $prev_beg0 line ".__LINE__."\n";
			print "$prog_tag GENE dont les limites de début ne sont pas dans un ordre croissant\n";
		    }
		    ($ID_0 eq '')or $ID_0 .= $separateur;
		    $ID_0 .= $gene[$cpt_deb_CDS];
		    push @id_numeric_gene_sens0, $cpt_deb_CDS;
		    print "$prog_tag OUI pour $ID_0\n";
		    if(($cpt_deb_CDS == $#{$beg[ $H_ind_T{'gene'} ]}) or ($cpt_deb_CDS == $dernier_deb_CDS_sens0))
		    {
			$cpt_seq++;
			print FSOR300_0SEQ ">contigues_up_seq_sens0_nbr_$cpt_seq dist_deb_trad ".(-$dist_fin_up_seq_deb_trad).": $ID_0".":\n".&print50($fasta,substr($sequence,$prev_beg0,$prev_end0 - $prev_beg0 + 1));
			print "$prog_tag on ecrit deb $prev_beg0, longueur ".($prev_end0 - $prev_beg0 + 1)." qui correspond a $prev_end0 - $prev_beg0 + 1 ligne ".__LINE__." gene : $ID_0\n";

			if($edit_limits_of_overlapped_seq_sens01){ $chaine_sens0 .= $prev_beg0.",".$prev_end0.";";}		
			if($edit_annotations_of_limits_overlapped){ print FSORANNOT_OVERLAPPED_SENS0 &get_string_annot($ID_0,\@id_numeric_gene_sens0);}
		    }
		    # exit;
		}
		else
		{
		    $cpt_seq++;
		    print FSOR300_0SEQ ">contigues_up_seq_sens0_nbr_$cpt_seq dist_deb_trad ".(-$dist_fin_up_seq_deb_trad).": $ID_0".":\n".&print50($fasta,substr($sequence,$prev_beg0,$prev_end0 - $prev_beg0 + 1));
		    print "$prog_tag on ecrit deb $prev_beg0, longueur ".($prev_end0 - $prev_beg0 + 1)." qui correspond a $prev_end0 - $prev_beg0 + 1 ligne ".__LINE__." gene : $ID_0\n";

		    if($edit_limits_of_overlapped_seq_sens01){ $chaine_sens0 .= $prev_beg0.",".$prev_end0.";"; }
		    
		    if($edit_annotations_of_limits_overlapped){ print FSORANNOT_OVERLAPPED_SENS0 &get_string_annot($ID_0,\@id_numeric_gene_sens0); }

		    $ID_0 = '';
		    @id_numeric_gene_sens0 = ();
		    $ID_0 .= $gene[$cpt_deb_CDS];
		    push @id_numeric_gene_sens0, $cpt_deb_CDS;
		    $prev_end0 = $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - 2 - $dist_fin_up_seq_deb_trad;
		    $prev_beg0 = $beg[ $H_ind_T{'gene'} ][$cpt_deb_CDS] - $dist_deb_up_seq_deb_trad - 1;
		    ($prev_beg0 >= 0)or $prev_beg0 = 0;	
		    
		    # si c'est le dernier gene
		    if(($cpt_deb_CDS == $#{$beg[ $H_ind_T{'gene'} ]}) or ($cpt_deb_CDS == $dernier_deb_CDS_sens0))
		    {
			$cpt_seq++;
			print FSOR300_0SEQ ">contigues_up_seq_sens0_nbr_$cpt_seq dist_deb_trad ".(-$dist_fin_up_seq_deb_trad).": $ID_0".":\n".&print50($fasta,substr($sequence,$prev_beg0,$prev_end0 - $prev_beg0 + 1));
			print "$prog_tag on ecrit deb $prev_beg0, longueur ".($prev_end0 - $prev_beg0 + 1)." qui correspond a $prev_end0 - $prev_beg0 + 1 ligne ".__LINE__." gene : $ID_0\n";
			if($edit_limits_of_overlapped_seq_sens01){ $chaine_sens0 .= $prev_beg0.",".$prev_end0.";"; }
		        if($edit_annotations_of_limits_overlapped){ print FSORANNOT_OVERLAPPED_SENS0 &get_string_annot($ID_0,\@id_numeric_gene_sens0); }
		    }		    		    
		}
		 
	    }
	
    }
    
    if($edit_limits_of_overlapped_seq_sens01) 
    {
	print FSORLIMSEQ "$chaine_sens0\n$chaine_sens1\n";
	close FSORLIMSEQ;
	print "$prog_tag FSORLIMSEQ file closed\n"; 
    }

    if($edit_annotations_of_limits_overlapped){


#      if($#ARGV < 2){
      if(not $bool_distinctfiles_sens0_1){
	system("cat $tmp_file_name[ 3 ] > $tmp_file_name[ 6 ]")  == 0 or die "$prog_tag [Error] system de 'cat $tmp_file_name[ 3 ] > $tmp_file_name[ 6 ]' failed line ".__LINE__.":$?\n";
	system("cat $tmp_file_name[ 4 ] >> $tmp_file_name[ 6 ]") == 0 or die "$prog_tag [Error] system de 'cat $tmp_file_name[ 4 ] > $tmp_file_name[ 6 ]' failed line ".__LINE__.":$?\n";
      }
      $bool_distinctfiles_sens0_1 and unlink ($tmp_file_name[ 3 ], $tmp_file_name[ 4 ]);

      #print FSORANNOT_OVERLAPPED_SENS1 "\n";
      close FSORANNOT_OVERLAPPED_SENS1;
      #print FSORANNOT_OVERLAPPED_SENS0 "\n";
      close FSORANNOT_OVERLAPPED_SENS0;

      if(not $not_verbose){ 
	foreach(@tmp_file_name[$#tmp_file_name - 2..$#tmp_file_name]){ print "$prog_tag $_ file closed\n"; } 
      }
    }

   # if($#ARGV < 2){
    if(not $bool_distinctfiles_sens0_1){
      system("cat $tmp_file_name[ 0 ] > $tmp_file_name[ 5 ]")  == 0 or die "$prog_tag [Error] system de 'cat $tmp_file_name[ 0 ] > $tmp_file_name[ 5 ]' failed line ".__LINE__.":$?\n";
      system("cat $tmp_file_name[ 1 ] >> $tmp_file_name[ 5 ]") == 0 or die "$prog_tag [Error] system de 'cat $tmp_file_name[ 1 ] > $tmp_file_name[ 5 ]' failed line ".__LINE__.":$?\n";
    }
    # to delete overlapped_seq_sens0 , overlapped_seq_sens1 
    $bool_distinctfiles_sens0_1 and unlink ($tmp_file_name[ 0 ], $tmp_file_name[ 1 ]);

    #print FSOR300_0SEQ "\n";
    close FSOR300_0SEQ;
    #print FSOR300_1SEQ "\n";
    close FSOR300_1SEQ;
    
    foreach(@tmp_file_name){   print "$prog_tag $_ file created\n"; }
 ##??? is this the end of the program??

    $cpt_seq = 0;
}

if($bool_premier_passage){ $bool_premier_passage = 0; }
}



sub get_limits($)
{
    my ($indice) = @_;
    #print "INDICE : $indice\n";
    (defined $indice) or die "$prog_tag [Error] Indice pas definit line ".__LINE__."\n";
    if(not defined $beg[ $H_ind_T{'gene'} ][$indice] )
    {
	die "$prog_tag [Error] Probleme no limit for $gene[$indice] line ".__LINE__."\n";
    }
    my $gene_beg = $beg[ $H_ind_T{'gene'} ][$indice];
    my $gene_end = $end[ $H_ind_T{'gene'} ][$indice];

    # print "gene_beg, gene_end $gene_beg, $gene_end\n";
    return ($gene_beg, $gene_end);
}


sub get_string_annot($\@)
{
   my ($id,$r_id_numeric_gene) = @_;
   my $chaine_annot = ">for_annot_$id\n";
   
   for my $numg(@$r_id_numeric_gene)
   {
       $chaine_annot .= "$gene[$numg] ";
       $chaine_annot .= ", prot_name $prot_name[$numg] : \n";

       if(not defined $prot_name[$numg] ){ die "$prog_tag [Error] prot_name non defni pour le gene $gene[$numg]\n";}
 #       if(defined $imp_linked_to[$numg])
#        {
#        	   for my $nri(0..$#{ $imp_linked_to[$numg] })
# 	   {
# 	       $chaine_annot .= "$imp_type_kind[ $imp_type[ $imp_linked_to[$numg][$nri][1] ]] , imp_sens $sens[ $H_ind_T{'imp'} ][ $imp_linked_to[$numg][$nri][1] ][ $imp_linked_to[$numg][$nri][0] ], $beg[ $H_ind_T{'imp'} ][ $imp_linked_to[$numg][$nri][1] ][ $imp_linked_to[$numg][$nri][0] ]..$end[ $H_ind_T{'imp'} ][ $imp_linked_to[$numg][$nri][1] ][ $imp_linked_to[$numg][$nri][0] ]";
# 	       defined($comment[ $H_ind_T{'imp'} ][$nri] ) and do { 
# 		   foreach(@{ $comment[ $H_ind_T{'imp'} ][$nri] }){
# 		       $chaine_annot .= ", comment_imp: $_\n"; 
# 		   }
# 	       }	      
# 	   }
#        }
#        else{ print "imp_linked_to of $gene[$numg]  NOT DEF\n"; }

      
       $chaine_annot .= "gene_sens $sens[ $H_ind_T{'gene'} ][$numg], "; 

       $chaine_annot .= "$beg[ $H_ind_T{'gene'} ][$numg]..$end[ $H_ind_T{'gene'} ][$numg], ";

        (defined $Nr_acc_prot[$numg]) and $chaine_annot .= "Nr_acc_prot $Nr_acc_prot[$numg], ";
#	(defined $seqdesc_title[ $H_ind_T{'prot'} ][$numg])and $chaine_annot .= "seqdesc_title $seqdesc_title[ $H_ind_T{'prot'} ][$numg], ";
	
#         (defined $beg[  $H_ind_T{'prot'} ][$numg][0]) and do {	
# 	    for my $nb_beg(0..$#{ $beg[  $H_ind_T{'prot'} ][$numg] }){
# 		$chaine_annot .= "$beg[ $H_ind_T{'prot'} ][$numg]..$end[ $H_ind_T{'prot'} ][$numg], ";
# 		(defined $sens[ $H_ind_T{'prot'} ][$numg])and $chaine_annot .= "prot_sens $sens[ $H_ind_T{'prot'} ][$numg], ";
# 	    }
# 	};

	defined($comment[ $H_ind_T{'gene'} ][$numg]) and $chaine_annot .= "\ngene comment :$comment[ $H_ind_T{'gene'} ][$numg]; ";
	
# 	defined($comment[ $H_ind_T{'prot'} ][$numg]) and do {
# 	    foreach(@{ $comment[ $H_ind_T{'prot'} ][$numg] }){
# 		$chaine_annot .= "\nprot comment :$_; ";
# 	    }
# 	};

#         if(defined $comment_cdregion[$numg])
# 	{
# 	    $chaine_annot .= "\ncomment_cdregion : $comment_cdregion[$numg]";
# 	}

#        	if(defined $region_linked_to[$numg])
# 	{
	   
# 	    $chaine_annot .= "\nRegion(s) :";	    
# 	    for my $nri(0..$#{ $region_linked_to[$numg] })	   
# 	    {	
# 		$chaine_annot .= "$region[$region_linked_to[$numg][$nri]],";					
# 	    }
	   
# 	}

        $chaine_annot .= "\n";
   }
   return ($chaine_annot);
			
}


# ************************************************************************************
# END WRITING IN FILES
# ************************************************************************************

exit;
