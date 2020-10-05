#!/usr/bin/perl -w

use Cwd;
use File::Find;
#use Number::Format qw(:subs :vars);
#$THOUSANDS_SEP   = '';
#$DECIMAL_POINT   = ',';
#$INT_CURR_SYMBOL = 'DEM';


# my $formatted = format_number($number);
  
# the following is COMMAND LINE how to call the program with the reqd arguments

# Ftest_rech_motif_fic_diff_ds_seq_promot_averees_vari_spacer_triTwoBoxes64multi.pl sma sco user_data/demo/SIGffRid_orthologs user_data/demo/sortie_sco_rmes.gaussien_8_3 user_data/demo/markov_mod_order7_SAV_complete_genome_NC_003155.txt user_data/demo/markov_mod_order7_SCO_complete_genome_NC_003888.txt 14 20 1 1 0 user_data/demo www wwnnw 0


# SIGffRid 0.01 : APRES article JOBIM

# pack useless?
# ref number instead of trinuc useless?

# VERSION SANS UTILISATION DE GRAPPE
# use Bonferoni criteria for select interesting words in R'MES results
# use less memory when lots of seq, more if few sequences

# 41: split seq set according to spacer in sub SORT...
#     remove if not enough seq 

# 43: can use for statistic evaluation of final motifs either upstream seq overrepresentation (as in previous versions) or promoter seq overrepresentation (intergenic upstream seq...)

# 44: we split set of sequences with very variable spacers in subset with only contigus authorized variation in sub SORT_IND (to avoid ggaatw8-30gtt for example)

# 45: we do not use uppercase letter to get spacers between trinucleotides to have consensus motif but only pos_tri array
#     we put uppercase letters for original rmes words shared by the two upstream orthologous seq
#
# 46: ADAPT to USE of GAP TRINUC

# 47: use of global var instead of parameters (if possible)

# 48: GAP trinuc using lib to generate trinuc

# 49: effective use of gapped trinuc

# 51: try of using boxes with 2 nucleotids in and only if it is associated to another (sum of letters >= 6)

# 52: minimal number of letter extended for evaluation was used considering letter extended; IN 52, we consider letters in motif given by MATRICE subprogr

# 53: we treat in a same sort call all the sets of sequences corresponding to a probability of extension for a letter lower than 0.2 at the position corresponding to the minimal proba; display of a warning if minimal proab used is higher than 0.2
#     directories are dedicated to motif files (results of the computing of ratio upstr/whole genome)
 
# 54: split for final motifs OK

# 55: deals motifs through directories (more clean) and some bugs fixed

# 56: explore every low probabilities for extensions

# 57: version used to see proba for a sigma factor (max reached when we select the best)

# 58: get annotations related to an INTERESTING motif and create the corresponding file in a dedicated directory
#     for the first read of fasta global annotation, create an index file to reach directly good annotation in file
#     by using SeqUtilities.pm personal library
#     SO: there are 2 directories: one for global annotation, the other for its indexation (byte)
#  good results, BAD proba 160106, corrected 180106

# 59?: use files to store and reload hashes

# 61: matrice_create2: compute regexp to take into account only letters with low proba
# sub matrice_create2 corrected
# correction of computed number of upstream seq (we must not take into account seq without nucleotids...) (nb_upstr_seq)
# and we did not take into account last seq, corrected...
# FANNOT corrected in sub rapport_criteria

# 63: sort of sorted_ind CORRECTED 
# use of MMO for LTR test
 
# 64: possible use of upstream fusionned seq and annot with several files (each), so seq and
# annot file must share syntax, and they are in same directory (different from prev version)

# 64multi: real treatment of every low proba (other versions badly thought, did not work)

# WE HAVE TO USE vec TO OPTIMIZE MEMORY FOR SPACERS STORAGE

# require 5.004;
# external LIB ********************************************************************************
use strict;
use Env qw(SIGFFRID_DIR);
# use Storable;              # to store data struct into a file for treatment of every low proba
                           # for extension in sub SORT
# use Math::MatrixReal;      # used because of exponentiation overload for matrix object
# use diagnostics;

# personal LIB ********************************************************************************
use lib $SIGFFRID_DIR.'libperl'; # we use following lib:
use GappedTrinucGenerator;       # qw(proc_gapped_trinuc_generator longueur);
use SysUtilities;                # qw(date tri_rapide clone);
use SeqUtilities;                # qw(complement_seq print50 create_index_file load_index_f);
# END LIB *************************************************************************************

$|++; # to display print in the executing order
# var for debug

my $prog_tag = '[Ftest_rech_motif_diff_ds_seq_promot_averees_vari_spacer_triTwoBoxes64multi.pl]';
my $time_long_or_not;
my $prev_time_long_or_not;
my $author_mail = 'fabrice.touzain\@anses.fr';

if($#ARGV < 6)
{
    die join("\n", "\nTo use this programme:",
	     "$0\n",
	     "\t<ID first species (involved in result file naming)>",
	     "\t<ID second species (involved in result file naming)(species of interest: real position used) >",
	     "\t<directory where are located files of ortholog sequences>",
	     "\t<formated rmes ouptut file>",
	     "\t<Markov_order analysis files for first species>",
	     "\t<Markov_order analysis files for second species>",
	     "\t<diff_min>",
	     "\t<diff_max>",
	     "\t<spacer :max difference allowed for a site in the given bacteria >",
	     "\t<spacer_shift : variation allowed for the same site in the 2 given ortholog seq>",
	     "\t<working_directory>",
	     "\t<seed1>",
	     "\t<seed2>",
	     "\t<boolean to tell if www special seed has already been treated>",
	     "\t<[file_search_string ]:> for the moment: not_fasta_",
    	     "\t<[file_search_string ]:> for the moment: -350_0_ \n"
	     
	);
}

# we are going to add 4 more arguments to be send to the script
# we will add : 1.  $diff_min: spacer min for 1 binding site (14) ARGV[6]
# 2.  $diff_max: spacer max for 1 binding site (20) ARGV[7]
# 3.  $spacer_shift_fic: max différence allowed for a same site in 1 bacteria,(1) ARGV[8]
# 4.  $spacer_shift: spacer variation allowed for the same site in 2 ortholog seq ([-spacer_shift..+spacer_shift] (1)  4 max ARGV[9]

## this is just to chk for parameters sent
my $numArgs = $#ARGV + 1;
print "\n \n Thanks, you gave me $numArgs arguments.\n";

 foreach my $argnum (0 .. $#ARGV) {
 
    print "$ARGV[$argnum]\n";
 
}

#system("echo \"\" | mail -s \"prog finished\" $author_mail");
#exit;

## end of chk

# to compute time of execution
my $debut_sigffrid = time();

# PARAMETERS FOR INPUT SEARCHED WORDS ************************************************************************
# boolean to use only words given with a positive score by RMES or not
my $positive_words_only = 0;
my $reg_zscore_searched;

if($positive_words_only)
{
    $reg_zscore_searched = '(\w+)\s+([\d\.inf]+)';
}
else
{
    $reg_zscore_searched = '(\w+)\s+([\d\.\-inf]+)';
}
# determine les mots dont l'ordre nous interesse (le fichier de diff MM est un fichier par valeur ordonnee dans chaque ordre, et les ordres sont DECROISSANTS)
my $ordre_mini_moins1 = 1; # minimal size of used words is equivalent to this value plus 2
my $ordre_maxi        = 6; # maximal size of used words is equivalent to this value plus 1

# boolean to know if we treat every combination of possible extensions with low proba (1) or only this with the lowest proba (0)
my $bool_use_every_low_proba = 0;

my %alpha = (0.001 => 0, 0.005 => 1, 0.01 => 2, 0.05 => 3, 0.1 => 4);

my $inf_zscore = 99999999; 
# ************************************************************************
# PARAMETERS FOR SELECTION OF MOTIFS ACCORDING TO LIKELYHOOD RATIO OF POISSON LOW
my $alpha_likelyhood_ratio = 0.05;
#                        seuil 0.001, 0.005, 0.01, 0.05, 0.1
my @POISSON_THRESHOLD_KHI2 = ( 10.83, 7.88 , 6.63, 3.84, 2.71 );
my $POISSON_THRESHOLD = $POISSON_THRESHOLD_KHI2[ $alpha{$alpha_likelyhood_ratio} ];
@POISSON_THRESHOLD_KHI2 = ();

# ************************************************************************
# PARAMETERS FOR SELECTION OF RMES WORDS, OUTPUT FILES AND DEBUG SELECTION
my $alpha_val = 0.005;
my @RMES_THRESHOLDS = ();
# score suit N(0,1)
# P(N(0,1)>t) = alpha/4^h ou h est la longueur du mot

# CORRECTED 200606

# threshold only for positive scores
# le seuil t aura pour valeur une de celles du tableau ci-dessous, en fonction des cas
# mots de longueur         3,    4,    5,    6,    7,    8       alpha
@{$RMES_THRESHOLDS[0]} = (4.16, 4.47, 4.75, 5.03, 5.29, 5.53); # seuil 0.001
@{$RMES_THRESHOLDS[1]} = (3.78, 4.11, 4.42, 4.71, 4.98, 5.24); # seuil 0.005
@{$RMES_THRESHOLDS[2]} = (3.60, 3.95, 4.27, 4.56, 4.85, 5.12); # seuil 0.01
@{$RMES_THRESHOLDS[3]} = (3.16, 3.54, 3.89, 4.22, 4.52, 4.80); # seuil 0.05
@{$RMES_THRESHOLDS[4]} = (2.95, 3.35, 3.72, 4.06, 4.37, 4.66); # seuil 0.1

# ************************************************************************
# ************************************************************************

# PARAMETERS FOR DISPLAY, OUTPUT FILES AND DEBUG SELECTION ***************
my $not_verbose         = 0; 
my $bool_print_cross    = 0;
my $bool_print_sp_align = 1;
my $bool_edit_fic_paire_ortho = 0;
# ************************************************************************

# PARAMETERS FOR MOTIF SELECTION ************************************

# boolean to know which kind of upstream seq will be used to compute scores (where regular expression motifs will be counted)
my $bool_intergenic_for_evaluation = 0; # 0 = we want to use only the upstream regions (upto -350) and 1= all the intergenic regions
my $motif_threshold            = 0.35; # 0.25; # taken into account only if FIXED TRESH
my $bool_fixed_motif_threshold = 0;
my $multipl_factor_motif_tresh = 3;    # interesting motif must be at less multipl_factor_motif_tresh
                                       # more present in upstream (or intergenic upstream) regions
# reading what the user gives as his choice
$bool_intergenic_for_evaluation = $ARGV[10];
$bool_intergenic_for_evaluation and $multipl_factor_motif_tresh = 5; # criteria must be less stringent for intergenic upstream because we have few 									sequences

# use for file naming according to threshold
$bool_fixed_motif_threshold or $motif_threshold = 'multi'.$multipl_factor_motif_tresh;

# minimal size for total number of letters in the two boxes, used to limit results when we use words of size 2 or 3. All the results with two boxes of only 6 letters will be removed. (associated with the following boolean)
my $mini_number_of_letters_in2trinuc   = 6; # 5 
my $maxi_number_of_letters_in2trinuc   = 6; 

# minimal size for total number of letters in -35 and -10 boxes
my $mini_number_of_letters_in2words   = 6;
my $bool_use_filter_for_mini_nb_of_comm_letters = 1; # to use the two previous limits (variables)

my $diff_min = $ARGV[6]; # default = 14  #12; # 14; 8;
my $diff_max = $ARGV[7];  # default = 20  #18; # 20; 14;
my $nb_seq_min_owning_a_word_in_one_ortholog_file = 2; # INUTILE??

# if we want to consider percentage (in the other case, we only fix the minimal number of occurrences)
my $bool_pc_actif = 0;       
# used only if bool_pc_actif == 0
my $nb_mini_occ   = 8; # 4       
# used only if bool_pc_actif == 1
# threshold to consider a motif interesting (20 sequences on 100 must have the motif) used only if bool_pc_actif == 1
my $percent_threshold = 1.5; 
# this value will be multiplied by the number of analysed files (number of couples of orthologs)
# used to limit number of non interesting results and speed treatment sometimes (use percentage or fixed thresold depending of bool_pc_actif)
my $min_nb_of_waited_sigma_motifs = $percent_threshold / 100; 
my $max_nb_of_waited_sigma_motifs = 6000;

# a result is consider as a good one if (with other criteria) total number of fixed letters for a regexp motif is lower than this filtre_nb_l_added_matrice_maxi (higher give bad statistics because sequences are in two few number) and higher than filtre_nb_l_added_matrice_mini (lower make us taking into account trinuc dyad: we loose only time)
my $filtre_nb_l_added_matrice_maxi = 14; 
my $filtre_nb_l_added_matrice_mini = 6;
# ICI: a ajouter pour diminuer complexite
# my $filtre_nb_l_in_1box            = 8; 

my $bool_degenerescence_max_motif = 0;

# pour la gestion des espacements variables entre boites
# shift autorise pour le comptage (et l'affichage) des resultats communs de doublets de trinucleotides)
my $spacer_shift_fic = $ARGV[8]; # default = 1; 
my @spacer_shift_fic = ();
# shift autorise dans l'espacement entre deux boites dans des sequences orthologues
my $spacer_shift = $ARGV[9]; # default = 1;  
my @spacer_shift = (0);

# ************************************************************************

# ************************************************************************
# PATH TO SET (DATA FOR STAT CRITERIA) 
# ************************************************************************
# we save in this array the two MM models file passed to the program, one for each bacteria
my @MM07_not_fasta_whole_seq_total = ($ARGV[4],$ARGV[5]) ;

# we get begin location to localize other ways according to this one
# ?? MUST BE DONE, CURRENT_DIR CAN BE DIFFERENT FROM prog_location (not the case for the moment)


my $working_dir = $ARGV[11];
my @file_search_string ;# this variable will store the prefix of the output files we get from extrait-seq-....pl script
			## for the moment this variable is " $file_search_string[0]= not_fasta_ && $file_search_string[1]= -350_0_ "
if (!defined($ARGV[15])) 
{ 
    @file_search_string = ('not_fasta_','-350_0_');
    print "$prog_tag file_search_string args not defined, set to ".join(',', @file_search_string).", line ".__LINE__."\n";

} # this is optional if not defined used default values
else
{
    @file_search_string= @ARGV[15..$#ARGV];
}

my $results_dir = $working_dir; ## right now its the same
# chomp($results_dir);
$results_dir .= "/SIGffRid_results/"; # the result directory where we want the results to be
if (! -d $results_dir) {mkdir $results_dir; print "$prog_tag Creation of $results_dir directory\n"}
# close CH;

# allow to accept both relatives and hard way to reach a directory
# my $results_dir_conditionnel = $results_dir;
# die "$results_dir toto\n";

my $way_to_upstr_seq_dir = ''; # path to reach upstream seq file(s) (or promoter seq file(s)) this is in working directory
my $deb_file = ''; # syntax of the beginnning of upstream seq files

my $way_to_annot_dir;

# if we have token seq from -350 to +10, shift_between_end_seq_and_deb_trad will be 10 (used to compute motif position according to traduction start)
my $shift_between_end_seq_and_deb_trad = 0; 

my $MOTIF_DIR;
my $way_to_annot_motif_dir;

# store tmp data (seq, annot, tmp files)
my $DATA_DIR;

if($bool_intergenic_for_evaluation) # ==1 ( i.e we chose only intergenic region)	
{
        $way_to_upstr_seq_dir = $working_dir."/";
   	$way_to_annot_dir = $working_dir."/"; 
    	$deb_file = $file_search_string[0]."promot_";
   	
	$MOTIF_DIR = $results_dir.'MOTIF_tresh'.$motif_threshold.'_';
   	$way_to_annot_motif_dir = $working_dir.'/MOTIF_ANNOT/';
    	$DATA_DIR = $results_dir.'DATA_DIR/';
}

else{ #default option: $bool_intergenic_for_evaluation==0 ( i.e we choose upstream region upto 350 bp )
	$way_to_upstr_seq_dir = $working_dir."/";
        # $file_search_string[0] = not_fasta_  and $file_search_string[1] = -350_0_
	if($file_search_string[1] !~ /_$/)
	{
	    $file_search_string[1] .= '_';
	}	
   	$file_search_string[1] =~/_([-\d]+)_/; 
	$shift_between_end_seq_and_deb_trad = $1;
       # these files are : not_fasta_-350_0_annot_overlapped_seq_sens1/0_bacid1.txt
	$way_to_annot_dir = $working_dir."/".$file_search_string[0].$file_search_string[1]."annot_overlapped_seq_sens01_";
        $deb_file = $file_search_string[0].$file_search_string[1];
	$deb_file .= "overlapped_seq_sens01_";
#	$way_to_annot_dir = $working_dir."/".$file_search_string[0].$file_search_string[1]."annot_overlapped_seq_sens01_";
#        $deb_file = $file_search_string[0].$file_search_string[1]."overlapped_seq_sens01_";

	$MOTIF_DIR = $results_dir.'MOTIF_UPSTR_tresh'.$motif_threshold.'_';
    	$way_to_annot_motif_dir = $working_dir.'/MOTIF_ANNOT_UPSTR/';
    	$DATA_DIR = $results_dir.'DATA_DIR_UPSTR/';
    	
    
}
# creation of directory of motifs file (to avoid to compute sevreral times R and LRT)
my $MOTIF_DIR_bact_min_nb_of_waited_sigma_motifs = '';

(-d $way_to_annot_motif_dir)or mkdir($way_to_annot_motif_dir);

if(! -d $DATA_DIR){
    mkdir($DATA_DIR);
}


my $way_to_whole_genome  =  $working_dir."/";
my $deb_file_whole_genome = $file_search_string[0]."whole_seq0_";

# file used for motif extension based on MM4

# my $Bernouilli_not_fasta_whole_seq_total = '../donnees_brutes/010205/bernouilli010205/bernouilli7_not_fasta_whole_seq0_SCO_total_date291004.txt';

# ************************************************************************
# PATH TO SET (DATA FOR STAT CRITERIA) 
# ************************************************************************


# CONSTANT used to know ranges analysed in sequence (used by GET2SEQFASTAMINMAX subprog)
my $LEFT_RIGHT_BOARDS_WIDTHS = 10;

# *****************PROBA
# parameter for proba computing
# my $seuil = 0.01;
# my $delta2_on_N = - log($seuil) / 2; # has to be multiplied by N, and then squarerooted so as to obtain real significative threshold

# PARAMETER GIVEN TO THE PROG *************************************

my $fic_mots_RMES = $ARGV[3]; # fichier resultant de la difference de fichiers de MM (les mots interessants sont positifs, tries dans l'ordre alphabetique, et chaque ordre est precede de la ligne # ordre n
# if($fic_diff =~ /^\//)
# {
  #  $results_dir_conditionnel = '';
# }
# else
# {
  #  $results_dir_conditionnel = $results_dir;
# }
(-e $fic_mots_RMES) or die "$prog_tag [Error] $fic_mots_RMES file does not exist, line ".__LINE__."\n";

open(FIC_DIFF,"< $fic_mots_RMES")or die "$prog_tag [Error] Impossible to open $fic_mots_RMES file:$!, line ".__LINE__."\n";

# ************************************************************************
# variables for names of results files
# ************************************************************************
my $entete_fic_res_1_pair = 'motifs_';
my @name_fic_sp = @ARGV[0..1]; # file to record align sequences of one of the orthologues so as to group related sequences of a same species
if($name_fic_sp[0] eq $name_fic_sp[1]){ $name_fic_sp[0] .= '1'; $name_fic_sp[1] .= '2'; }


# ************************************************************************
# VAR for directory and files used by Storable
# ************************************************************************
# my $suffix_store_file = '.dat';
# my $cpt_store_file    = -1; # count number of stored files to avoid to erase file of another "stack" of the same couple of trinuc
# my $nb_file_to_treat  = 0;
# print "cpt_store_file initisalise a $cpt_store_file\n";

# ************************************************************************

if($bool_pc_actif){  $MOTIF_DIR .= 'pc_';                           }
else              {  $min_nb_of_waited_sigma_motifs = $nb_mini_occ; }

for my $s(@name_fic_sp)
{
    my $complete_dir_name = $MOTIF_DIR.$s.'_'.$min_nb_of_waited_sigma_motifs;
    if(not -e $complete_dir_name)
    {
	mkdir($complete_dir_name)or die "$prog_tag [Error] Impossible to create $complete_dir_name:$!, line ".__LINE__."\n";
    }
    print "$prog_tag $complete_dir_name directory created, line ".__LINE__."\n";
}

# END PARAMETER GIVEN TO THE PROG *********************************

# Hash table used to compute pseudo-probability to lead extension
my %H_Bernouilli_extension = (); # keys: trinuc, 
my %H_MM4_extension = (); # keys: trinuc, 
# first col: letter (before or after)
# second col: P(letter) before
# third col: P(letter) after

my $bool_not_BERNOUILLI = 1;

# to associate to a letter a table indice
my %H_l = ('a' => 0,
	   'c' => 1,
	   'g' => 2,
	   't' => 3,
	   'n' => 4,
	   'A' => 0,
	   'C' => 1,
	   'G' => 2,
	   'T' => 3,
	   'N' => 4);

# to get letter associated to a table indice
my @rev_H_l         = ('a','c','g','t','n'); 
my $bool_www_treated = $ARGV[14]; ## TO DOOO
# my @VAR_trinuc     = ('www3','wwww','wnnw','wnnnw','wnw','wwnnw','wnnww','wnwnw','wwnw','wnww'); # store every 'regexp' needed to generate every declination of 'trinuc'
my @VAR_trinuc    = @ARGV[12 .. 13]; # store every 'regexp' needed to generate every declination of trinuc 

######## these are the names of files
my $align_id1 = "align_".$diff_min."_".$diff_max."_".$spacer_shift_fic."_".$spacer_shift."_".$VAR_trinuc[0]."_".$VAR_trinuc[1]."_".$name_fic_sp[0].".txt";
my $align_id2 = "align_".$diff_min."_".$diff_max."_".$spacer_shift_fic."_".$spacer_shift."_".$VAR_trinuc[0]."_".$VAR_trinuc[1]."_".$name_fic_sp[1].".txt" ;



($#VAR_trinuc < 10)or die "$prog_tag [Error] You can not use more than 10 regexp (trinuc) to group interesting results, line ".__LINE__."\n";

my @table_trinuc   = ();      # store every declination of trinuc used for grouping
my @trinuc_lengths = ();      # store trinuc lengthes (first field) to associate regexp indices 
                              # (of VAR_trinuc) (second field = anonymous tables)

# foreach trinuc in table trinuc, give for its indice (one letter of the chain), the indice of corresponding regexp in VAR_trinuc (FIXE limit number of regexp at ten values)
my $str_trinuc_give_VAR = '';

# call GAPPED_TRINUC_GENERATOR in gapped_trinuc_generator library
# CAUTION: one table is sorted:  other table recorded after must not be recorded before otherelse, we will have bad relations between tables
&proc_gapped_trinuc_generator(\@VAR_trinuc,\@table_trinuc,\@trinuc_lengths,\$str_trinuc_give_VAR);


# for(my $i = 0; $i <= $#table_trinuc; $i++)
# {
#     my $tmp_substr = substr($str_trinuc_give_VAR, $i, 1);
#    ($not_verbose)or print "num $i, tri $table_trinuc[$i], re $tmp_substr\n";
# }
 # exit;

my @tmp_VAR_trinuc   = @VAR_trinuc; # store every 'regexp' corresponding to trinuc

# store nb of fixed letters in a trinuc (avoid to compute it n times)
my @nbl = (0) x @VAR_trinuc;
for(my $i = 0; $i <= $#VAR_trinuc; $i++){ 
    while($VAR_trinuc[$i] =~ /w/g){ $nbl[$i]++; }
    ($not_verbose)or print "$i VAR_trinuc $VAR_trinuc[$i], nbl $nbl[$i]\n";
}

# verif @nbl
# for(my $i = 0; $i <= $#nbl; $i++)
# {
  #  print "$VAR_trinuc[$i] $nbl[$i]\n";
# }
# exit;

# ******************************************************************************
# CREATION OF REGEXP USED TO HAVE TRINUC WITH GAPS TO GROUP MOTIFS
# ******************************************************************************
for my $re(@tmp_VAR_trinuc)
{
    # print "re $re ";
    &REMPLACE_W_PAR_ANTISLASHW(\$re);
    # print "devient $re\n";

    # print "re $re ";
    &REMPLACE_N_PAR_ANTISLASHW(\$re);
    # print "devient $re\n";
    
}
# END CREATION OF REGEXP USED TO HAVE TRINUC WITH GAPS TO GROUP MOTIFS
# ******************************************************************************

# ******************************************************************************
# DSIPLAYING FOR VERIFICATION
# ******************************************************************************
for(my $i = 0; $i <= $#tmp_VAR_trinuc; $i++)
{
    print "i $i, VAR_trinuc $VAR_trinuc[$i], tmp_VAR_trinuc $tmp_VAR_trinuc[$i]\n";
}
# exit;
# @VAR_trinuc   = (); # store every 'regexp' needed to generate every declination of trinuc

print "tab_trinuc has ".scalar(@table_trinuc)." entries\n";
print "tab_trinuc\n";
my $cptcell = 0;
foreach(@table_trinuc)
{
    print $_,' ';
    (++$cptcell % 4 == 0)and print "\n";
}
# ******************************************************************************

# ******************************************************************************
# 57 ****************************** to deal with files record for hash
# ******************************************************************************
# my $useDiskStorage = 1; # active BD use (can be deactivated if small bacteria (faster))
# my $storage_DIR    = $MOTIF_DIR;
# $storage_DIR       =~ s/MOTIF/STORAGE_DIR/;
# my $storage_DIR   .= join('_', @name_fic_sp, $min_nb_of_waited_sigma_motifs, @VAR_trinuc);
# (-d $storage_DIR)or mkdir($storage_DIR)or die "Can not create storage_DIR $storage_DIR\n";

# my $key2_balise = 'key2';

# 57 ****************************** end to deal with files record for hash *****
# ******************************************************************************

# my @table_trinuc = ('aaa','aac','aag','aat',
	#	    'aca','acc','acg','act',
	#	    'aga','agc','agg','agt',
	#	    'ata','atc','atg','att',
	#	    'caa','cac','cag','cat',
	#	    'cca','ccc','ccg','cct',
	#	    'cga','cgc','cgg','cgt',
	#	    'cta','ctc','ctg','ctt',
	#	    'gaa','gac','gag','gat',
	#	    'gca','gcc','gcg','gct',
	#	    'gga','ggc','ggg','ggt',
	#	    'gta','gtc','gtg','gtt',
	#	    'taa','tac','tag','tat',
	#	    'tca','tcc','tcg','tct',
	#	    'tga','tgc','tgg','tgt',
	#	    'tta','ttc','ttg','ttt'
	#	   );

# table to determine which trinuc couple to treat in the same time to limit memory usage according to GC percent ??
# my @tab_repartition;

# H table for trinuc H tables
my %H_trinuc = ();

# initialization if only trinuc words
my $count_trinuc = 0;
foreach(@table_trinuc)
{
    $H_trinuc{ $_ } = $count_trinuc++; # we associate key 'aaa' to indice in 
# $not_verbose or print $rev_H_l[$first],' ',$rev_H_l[$second],' ',$rev_H_l[$third]," traite\n";
}

#die "table_trinuc $table_trinuc[$H_trinuc{'aaa'}]\n";
# print "verif H_trinuc\n";
# foreach my $cle(sort keys %H_trinuc) 
# {
#    print "$cle = $H_trinuc{$cle}\n";
# }
# print "end verif H_trinuc\n";
# exit;

my %MM0 = ();

foreach my $sp(0..$#name_fic_sp)
{
    # !!!!! this file is treated as Bernouilli file (do not change results for MM but change results for Bernouilli!!)

    open(BERNOUILLI,"< $MM07_not_fasta_whole_seq_total[$sp]") or die "$prog_tag [Error] Impossible to open $MM07_not_fasta_whole_seq_total[$sp] :$!, line ".__LINE__."\n";
    
    $bool_not_BERNOUILLI = 1;
    my $treated_order;
    my $bool_MM0_treatment = 0;
    # Bernouilli H design according to Bernouilli file
    while(my $line = <BERNOUILLI>)
    {
	if($line =~ /^\#ordre 2/){ $bool_not_BERNOUILLI = 1; }
	elsif($line =~ /^\#ordre 0/){ $bool_MM0_treatment = 1; }

	# order 3 used for motif extension
	# order 1 used for motif proba for LRT test (likelyhood ratio test)
	elsif($line =~ /^\#ordre ([31])/){ $bool_not_BERNOUILLI = 0; $treated_order = $1; $bool_MM0_treatment = 0; }
	elsif($line =~ /^\#ordre 4/){ last; }
	elsif($bool_MM0_treatment){
	    $line =~ /^(\w)\s+([\d\.]+)/;
	    $MM0{$1}[$sp] = $2;
	}
	elsif( $bool_not_BERNOUILLI ){ next; }
	else{
	    $line =~ /^(\w+)\s+([\d\.]+)/;
	    # if(! defined $1){ die "line $line dol1 non def, line ".__LINE__."\n"; }
	    # proba for letter before trinuc
	    #               {  trinuc        }[ id bacterie ][ 0=before, 1=after ][ indice (letter)        ]
	    $H_Bernouilli_extension{ substr($1,1,$treated_order) }[ $sp ][ 0 ][ $H_l{ substr($1,0,1) } ] = $2; 
	    # proba for letter after trinuc
	    $H_Bernouilli_extension{ substr($1,0,$treated_order) }[ $sp ][ 1 ][ $H_l{ substr($1,$treated_order,1) } ] = $2; 
	    
	}
    }
    
    # verif %H_MM4_extension
    # we create MM4 hash from Bernouilli hash
    foreach my $trinuc( keys %H_Bernouilli_extension)
    {
	foreach my $b(0..$#{ $H_Bernouilli_extension{$trinuc} } )
	{
	    foreach my $pos(0..$#{ $H_Bernouilli_extension{$trinuc}[$b] } )
	    {
		
		my $total = 0;
		foreach my $l(0..$#{ $H_Bernouilli_extension{$trinuc}[$b][$pos] } )
		{
		    # print "$trinuc : lettre $l pos $pos val $H_Bernouilli_extension{$trinuc}[$pos][$l]\n";
		    
		    # proba for letter before trinuc
		    #               {  trinuc        }[ 0=before, 1=after ][ indice (letter)        ]
		    $total += $H_Bernouilli_extension{ $trinuc }[$b][ $pos ][ $l ];  
		}
		# print "total vaut $total\n";
		foreach my $l(0..$#{ $H_Bernouilli_extension{$trinuc}[$b][$pos] } )
		{
		    # proba for letter after trinuc
		    $H_MM4_extension{ $trinuc }[$b][ $pos ][ $l ] = $H_Bernouilli_extension{ $trinuc }[$b][ $pos ][ $l ] / $total; 
		}	    
	    }
	}
    }
    
    close BERNOUILLI;
    print "$prog_tag START VERIF %H_MM4_extension, line ".__LINE__.":\n";
    foreach my $trinuc( keys %H_Bernouilli_extension)
    {
	foreach my $pos(0..$#{ $H_Bernouilli_extension{$trinuc} } )
	{
	    foreach my $l(0..$#{ $H_Bernouilli_extension{$trinuc}[$pos] } )
	    {
		# Bernouilli
		print "trinuc $trinuc : lettre $l pos $pos val $H_Bernouilli_extension{$trinuc}[$pos][$l]\n";
		
		# MM4
		print "$trinuc : lettre $l pos $pos val $H_MM4_extension{$trinuc}[$pos][$l]\n";
	    }
        }
    }
    print "$prog_tag END VERIF %H_MM4_extension:\n";    
}

%H_Bernouilli_extension = ();


# *****************END PROBA

# ******************************************************************************
# *****************IMPRESSION INFOS
# ******************************************************************************

my $infos = "PARAMETRES: diff_min $diff_min diff_max $diff_max nb_mini_occ $nb_mini_occ spacer_shift_fic $spacer_shift_fic spacer_shift $spacer_shift, bool_intergenic_for_evaluation $bool_intergenic_for_evaluation, ordre_mini_moins1 $ordre_mini_moins1 ordre_maxi $ordre_maxi, bool_use_every_low_proba $bool_use_every_low_proba, mini_number_of_letters_in2trinuc $mini_number_of_letters_in2trinuc,  mini_number_of_letters_in2words $mini_number_of_letters_in2words, bool_use_filter_for_mini_nb_of_comm_letters $bool_use_filter_for_mini_nb_of_comm_letters, motif_threshold $motif_threshold, bool_fixed_motif_threshold $bool_fixed_motif_threshold, multipl_factor_motif_tresh $multipl_factor_motif_tresh, positive_words_only $positive_words_only, min_nb_of_waited_sigma_motifs $min_nb_of_waited_sigma_motifs, max_nb_of_waited_sigma_motifs $max_nb_of_waited_sigma_motifs\nVAR_trinuc @VAR_trinuc\nCALL: $0 @ARGV\n";

$infos .= "DATE: ".&date('/');
print $infos."\n";

my @tab_ID_fic = (); # associate ID (indice) to a file identifiant

# ******************************************************************************
# TABLE initialization for authorized shift in two upstream seq of orthologs
# (0 already recorded)
for my $variat(1..$spacer_shift)
{
    push @spacer_shift, -$variat;
    push @spacer_shift, +$variat;
}
# ******************************************************************************

# ******************************************************************************
# TABLE initialization for authorized shift in every upstream seq of a same b
for my $variat(reverse 1..$spacer_shift_fic)
{
    push @spacer_shift_fic, -$variat;
}
push @spacer_shift_fic, 0;
for my $variat(1..$spacer_shift_fic)
{
    push @spacer_shift_fic, +$variat;
}
# ******************************************************************************

my $directory_for_ortho_seq = $ARGV[2]; # dir where orthologus seq files are present
chomp($directory_for_ortho_seq);
(-d $directory_for_ortho_seq) or die "$prog_tag [Error] $directory_for_ortho_seq directory does not exist, line ".__LINE__."\n";
# we add slash at the end if not existing
if ($directory_for_ortho_seq!~/\/$/) {$directory_for_ortho_seq =~ s/([^\/]+)$/$1\//};

my %mots_surrep   = ();        # stocke les mots surrepresentes dans les regions intergeniques
                               # entre genes divergents par rapport aux regions intergeniques entre genes convergents
my @t_mots_surrep = ();        # idem mais ordonne

my $bool_ordre_interessant = 0;
my $length_RMES_w          = 0;

# ?? GO ON COMMENTING FROM HERE
# we get all the overrepresented words according to the file corresponding to the ratio between the probabilities 
# in two MM files (between divergent genes on between convergent genes)
my $tmp_zscore = undef;
while(<FIC_DIFF>)
{
    if(/^\#ordre $ordre_maxi/)
    { 
		$length_RMES_w = $ordre_maxi+1;
		$bool_ordre_interessant = 1; 
		($not_verbose)or print "RMES threshold $length_RMES_w letters $RMES_THRESHOLDS[$alpha{$alpha_val}][$length_RMES_w-3]\n";
    }
    elsif(/^\#ordre $ordre_mini_moins1/)
    {
		last;
    }
    elsif(/^\#ordre (\d+)/)
    {
		$length_RMES_w = $1+1;
		($not_verbose)or print "RMES threshold $length_RMES_w letters $RMES_THRESHOLDS[$alpha{$alpha_val}][$length_RMES_w-3]\n";
    }    
    elsif(! $bool_ordre_interessant)
    { 
		next; 
    }
    # ne concerne que les mots positifs puisqu'il n'y a pas de signe - avant le chiffre
    # elsif(/(\w+)\s+([\d\.]+)/)
    elsif(/$reg_zscore_searched/)
    {
		$tmp_zscore = $2;

		# case of infiny zscore
		foreach($tmp_zscore)
		{
			/^inf$/  and do{ $tmp_zscore =  $inf_zscore; next; };
			/^-inf$/ and do{ $tmp_zscore = -$inf_zscore; next; };
			/\d/ and next;
			die "$prog_tag [Error] $_ zscore not recognized as number or inf or -inf, line ".__LINE__."\n";
		}
		
		if(($tmp_zscore > $RMES_THRESHOLDS[$alpha{$alpha_val}][$length_RMES_w-3])or
		   ((! $positive_words_only)and($tmp_zscore < -$RMES_THRESHOLDS[$alpha{$alpha_val}][$length_RMES_w-3]))) 
		{
		    # record score resulting in division between sequences in different MM files
		    $mots_surrep{$1} = $tmp_zscore; 
		    push @t_mots_surrep, $1;
		}
    }
    # else{  print "Value less than 1 for word at line $_ in file, line ".__LINE__." in prog\n";   }
}
close(FIC_DIFF);

@RMES_THRESHOLDS = ();

# my @ind_lnUP_val_lnWHOL = (); # ??? store index: to each upstrem seq line (indice), we give a line in every upstream seq file

# verif t_mots_surrep
# print "t_mots_surrep @t_mots_surrep\n";
# exit;

# system("rm -f motif*.txt");


# my @name_ID_file = (); # only used to create ID files in order to avoid creating them again
# my @name_ID_file  = (); # table of table which record for each species (first indice) all the seq ID
# my @name_SEQ_file = (); # table which record seq file name where are searched regexp_motifs

# we do not compute ratio to ponderate ratio scores for motifs because when we have palyndromique, must we compute in the two sens for global count?
my @upstr_nt         = (0)x@name_fic_sp;
my @total_nt         = (0)x@name_fic_sp;
my @rapports         = (0)x@name_fic_sp;
my @rapports_display = (0)x@name_fic_sp;
my @nb_upstr_seq     = (0)x@name_fic_sp;
my @nb_whole_seq     = (0)x@name_fic_sp;

my @index_annot_f       = (); # names of annotation files foreach bacteria
my @annot_byte_position = (); # store 2 tables. Each one is a table of every byte positions of annotations
                              # in annotation files...
my %MM0_upstr           = (); # MM0 hash to associate a proba to a letter in upstr seq got
my ($handle_ref_glob00, $handle_ref_glob01, $handle_ref_glob10, $handle_ref_glob11) = ('')x4;
# global var to store SEQ and ID glob
my ($handle_ref_glob_sp0, $handle_ref_glob_sp1) = ('')x2; 

for my $sp(0..$#name_fic_sp){ # 1st for 

    # $name_ID_file[$sp]  = "ID_$name_fic_sp[$sp]".".txt";
    # $name_SEQ_file[$sp] = "SEQ_$name_fic_sp[$sp]".".txt";

    # print '\n prog_location '.$results_dir;
    # print 'ID_'.$name_fic_sp[$sp].".txt\n";

 
    # END creation of index files for annotation files (for faster searching and getting) ****************
    # my @regexp_lnWHOL_lnUP = ('','');

    if(not -e $DATA_DIR.'ID_'.$name_fic_sp[$sp].".txt")
    { # 1st if
	print "$prog_tag ".$DATA_DIR.'ID_'.$name_fic_sp[$sp].".txt file not found, line ".__LINE__."\n";
	my @file_all_promot = ();
		find sub { 
			if($File::Find::name=~ /$deb_file$name_fic_sp[$sp]/){
			  push @file_all_promot, $File::Find::name; 
			}
		}, "$way_to_upstr_seq_dir";
		
		print "Promoter file used:".join(',', @file_all_promot)."\n";
		print "They have to contain '$deb_file$name_fic_sp[$sp]' string and be in $way_to_upstr_seq_dir directory\n";
		
		# open(FICSEQALLPROMOT,'ls '.$way_to_upstr_seq_dir.$deb_file."$name_fic_sp[$sp]"."* |")or die "Impossible to open ls command line ".__LINE__.":$!\n";
		if($sp){

		  # store only SEQ promot each one on a single line
		  open(FSORSEQ1,  "+> $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt")or die "$prog_tag [Error] Impossible to create $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt line ".__LINE__.": $!\n";
		  # store only ID of promot, each one on a single line
		  open(FSORID1,   "+> $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt") or die "$prog_tag [Error] Impossible to create $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt line ".__LINE__.": $!\n";
		  $handle_ref_glob_sp0 = $handle_ref_glob10 = *FSORSEQ1; 
		  $handle_ref_glob_sp1 = $handle_ref_glob11 = *FSORID1;
		}
		else{
		  # store only SEQ promot each one on a single line
		  open(FSORSEQ0,  "+> $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt")or die "$prog_tag [Error] Impossible to create $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt line ".__LINE__.": $!\n";
		  # store only ID of promot, each one on a single line
		  open(FSORID0,   "+> $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt") or die "$prog_tag [Error] Impossible to create $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt line ".__LINE__.": $!\n";
		  $handle_ref_glob_sp0 = $handle_ref_glob00 = *FSORSEQ0; 
		  $handle_ref_glob_sp1 = $handle_ref_glob01 = *FSORID0;
		}
	print "$prog_tag $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt file created\n";
	print "$prog_tag $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt file created\n";		  
		
		# store only nb of nucleotides in all promot seq and in whole seq
		open(FSORNBNUC,"> $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]"."_nt.txt") or die "$prog_tag [Error] Impossible to create $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]"."_nt.txt line ".__LINE__.": $!\n";
	
		$MM0_upstr{a}[$sp] = $MM0_upstr{c}[$sp] = $MM0_upstr{g}[$sp] = $MM0_upstr{t}[$sp] = 0;
		# CAUTION: CAS DE PLUSIEURS FICHIERS D'ANNOT NON ENCOORE TRAITE ???
		foreach my $fic(@file_all_promot){ # while1
	
	
		    print "We use $fic file\n";
		    # next if($fic =~ /_total/);
		    open(FPROM,"< $fic")or die "$prog_tag [Error] Impossible to open $fic: $!, line ".__LINE__."\n";
		    # print FSOR "\n Fichier sens 0:\n";
		    # print "fic $fic traite\n";
		    # $fic =~ s/.*?([^\/]+)$/$1/;
		    
		    my $seq = '';
	
		   # my $regexp_lnWHOL_lnUP = '';
	
		    while(my $l = <FPROM>){ # while2
	
			if($l =~ /^>/){
				
			    # $regexp_lnWHOL_lnUP[$sp] .= "$1|";
			    ($seq ne '')and $nb_upstr_seq[$sp]++;
			    ($. == 1)or print $handle_ref_glob_sp0 "\n";
			    $seq = '';
			    print $handle_ref_glob_sp1 $l;
			}
			# elsif($l =~ /^>/){
			  #  die "promot $l not recognized line ".__LINE__."\n";
			# }
			else{
	
			    chomp($seq = $l);
			    print $handle_ref_glob_sp0 lc($seq);
			    $upstr_nt[$sp] += length($seq);
			    for(split //, $seq){ $MM0_upstr{$_}[$sp]++; }
			    # print "seq $seq ajoutee\n";
			}
	
		     } # end  of while2
	
		    # treatment for last seq (because there is no '>' after...)
		    ($seq ne '')and $nb_upstr_seq[$sp]++;
		    print $handle_ref_glob_sp0 lc($seq);
		    $upstr_nt[$sp] += length($seq);
		    for(split //, $seq){ $MM0_upstr{$_}[$sp]++; }
		    # $seq = '';
		  
	
		    # die "regexp_lnWHOL_lnUP $regexp_lnWHOL_lnUP\n";
	
		    # delete last '|' symbol
		    # chop $regexp_lnWHOL_lnUP[$sp];
		    # $regexp_lnWHOL_lnUP = $regexp_lnWHOL_lnUP; # ???
		    # print "avant GREP ".__LINE__."\n";
		    # system("grep -n '$regexp_lnWHOL_lnUP' \> GREP.txt")or die "PB sys line ".__LINE__."\n";
		    # open(GREP, "grep -n '$regexp_lnWHOL_lnUP' ??? |")or die "Can not call GREP line ".__LINE__."\n";
		    # print "apres GREP ".__LINE__."\n";
		    # ???
		    # open(LNUP_LNWHOL, "> $DATA_DIR".'lnUP_lnWOHL_ind_'."$name_fic_sp[$sp]".".txt")or die "Can not open $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]"."_nt.txt line ".__LINE__."\n"; 
		    # while(<GREP>){
			# /^(\d+):/ and do {
			  #  print LNUP_LNWHOL $1.' ';
			  #  push @{ $ind_lnUP_val_lnWHOL[$sp] }, $1; # ??? A MODIF
			  #  print "boucle ".__LINE__."\n";
			  #  next;
			# };
			# die "$_ unknown case for line $_, line ".__LINE__."\n";
		   #  }
		   #  print "boucle GREP passee ".__LINE__."\n";
		   #  close GREP;
	
		    close FPROM;
		} # end of while1

		$nb_upstr_seq[$sp] or die "$prog_tag [Error] No seq for promot or upstream seq, for sp $sp line ".__LINE__."\n";
		$MM0_upstr{a}[$sp] /= $upstr_nt[$sp];
		$MM0_upstr{c}[$sp] /= $upstr_nt[$sp];
		$MM0_upstr{g}[$sp] /= $upstr_nt[$sp];
		$MM0_upstr{t}[$sp] /= $upstr_nt[$sp];
	
		# close(FICSEQALLPROMOT)  or die "Impossible to close FICSEQALLPROMOT line ".__LINE__.":$!";
	
		# we search for total number of nucleotids
		# ***************************************************************
		
	
		$total_nt[$sp] = 0;

	print "$prog_tag We read $way_to_whole_genome".$deb_file_whole_genome.$name_fic_sp[$sp].".txt file, line ".__LINE__."\n";
	open(WHOLESEQ,'<', $way_to_whole_genome.$deb_file_whole_genome.$name_fic_sp[$sp].".txt")
	    or die "$prog_tag [Error] Impossible to open $way_to_whole_genome".$deb_file_whole_genome.$name_fic_sp[$sp].".txt: $!, line ".__LINE__."\n";
	
	while(my $seqw = <WHOLESEQ>){
	    
	    $total_nt[$sp] += length($seqw);
	    ($seqw ne '')and $nb_whole_seq[$sp]++;
	    # ($not_verbose)or print "length line ".length($seqw)."\n";
	} # end of while
	close WHOLESEQ;
		
 
		# ***************************************************************
	
		print FSORNBNUC "count upstr seq nt, total nt, nb_upstr_seq, nb_whole_seq: $upstr_nt[$sp] $total_nt[$sp] $nb_upstr_seq[$sp] $nb_whole_seq[$sp] $MM0_upstr{a}[$sp] $MM0_upstr{c}[$sp] $MM0_upstr{g}[$sp] $MM0_upstr{t}[$sp]\n";
		# close FSORSEQ.$sp;
		# close FSORID.$sp;
	close FSORNBNUC;
	print "$prog_tag $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]"."_nt.txt file created\n";

    }# 1st if ends here
    else
    {
	my $fsortnbnuc_f = $DATA_DIR.'SEQ_'."$name_fic_sp[$sp]"."_nt.txt";
	print "$prog_tag ".$DATA_DIR.'ID_'.$name_fic_sp[$sp].".txt file found, we open $fsortnbnuc_f file, line ".__LINE__."\n";

	my $first_line = undef;
	while(not defined $first_line)
	{
	    open(FSORNBNUC,'<',$fsortnbnuc_f)or die "$prog_tag [Error] Impossible to open $fsortnbnuc_f:$!, line ".__LINE__."\n";
	    
	    $first_line = <FSORNBNUC>;
	    if(not defined $first_line)
	    {
		warn "$prog_tag [Error] first line ($first_line) not defined? in $fsortnbnuc_f file, line ".__LINE__.", probably another Ftest_rech script is creating it, we try to open it again for reading after 2 minutes\n";
		sleep 7200;
	    }
	}
		
		if($first_line =~ /\: (\d+) (\d+) (\d+) (\d+) ([\d\.]+) ([\d\.]+) ([\d\.]+) ([\d\.]+)[\s\t]*$/)
		{
			$upstr_nt[$sp]     = $1;
			$total_nt[$sp]     = $2;
			$nb_upstr_seq[$sp] = $3;
			$nb_whole_seq[$sp] = $4;
			$MM0_upstr{a}[$sp] = $5;
			$MM0_upstr{c}[$sp] = $6;
			$MM0_upstr{g}[$sp] = $7;
			$MM0_upstr{t}[$sp] = $8;
			defined $1 or die "$prog_tag [Error] dol1 not defined, line ".__LINE__."\n";
			defined $2 or die "$prog_tag [Error] dol2 not defined, line ".__LINE__."\n";
			defined $3 or die "$prog_tag [Error] dol3 not defined, line ".__LINE__."\n";
			defined $4 or die "$prog_tag [Error] dol4 not defined, line ".__LINE__."\n";
			defined $5 or die "$prog_tag [Error] dol5 not defined, line ".__LINE__."\n";
			defined $6 or die "$prog_tag [Error] dol6 not defined, line ".__LINE__."\n";
			defined $7 or die "$prog_tag [Error] dol7 not defined, line ".__LINE__."\n";
			defined $8 or die "$prog_tag [Error] dol8 not defined, line ".__LINE__."\n";
		}
		else
		{
			die "$prog_tag [Error] first line ($first_line) not recognized by regexp in $fsortnbnuc_f file, line ".__LINE__."\n";
		}
		print "pour sp $sp, upstr_nt total_nt nb_upstr_seq nb_whole_seq MM0_upstr a MM0_upstr c MM0_upstr g MM0_upstr t: $1 $2 $3 $4 $5 $6 $7 $8 trouves ligne ".__LINE__."\n";
		close FSORNBNUC;
    } # else for 1st if ends

	if($sp){
	  
	  # store only SEQ promot each one on a single line
	  open(FSORSEQ1,  "< $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt")       or die "$prog_tag [Error] Impossible to read $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt line ".__LINE__.": $!\n";
	  # store only ID of promot, each one on a single line
	  open(FSORID1,   "< $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt")        or die "$prog_tag [Error] Impossible to read $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt line ".__LINE__.": $!\n";
	  $handle_ref_glob10 = *FSORSEQ1;
	  $handle_ref_glob11 = *FSORID1;
	}

	else{
	  # store only SEQ promot each one on a single line
	  open(FSORSEQ0,  "< $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt")       or die "$prog_tag [Error] Impossible to read $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]".".txt line ".__LINE__.": $!\n";
	  # store only ID of promot, each one on a single line
	  open(FSORID0,   "< $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt")        or die "$prog_tag [Error] Impossible to read $DATA_DIR".'ID_'."$name_fic_sp[$sp]".".txt line ".__LINE__.": $!\n";
	  $handle_ref_glob00 = *FSORSEQ0;
	  $handle_ref_glob01 = *FSORID0;
	}

#	open(FSORID.$sp,"< $DATA_DIR"."$name_ID_file") or die "Impossible to open $DATA_DIR"."$name_ID_file:$!";
#	@{ $name_ID_file[$sp] } = <FSORID.$sp>;
#	close FSORID.$sp;
 
	# open(LNUP_LNWHOL, "> $DATA_DIR".'lnUP_lnWOHL_ind_'."$name_fic_sp[$sp]".".txt")or die "Can not open $DATA_DIR".'SEQ_'."$name_fic_sp[$sp]"."_nt.txt line ".__LINE__."\n"; 
	# print "avant first_line ".__LINE__."\n";
	# my $first_line = <LNUP_LNWHOL>;
	# @{ $ind_lnUP_val_lnWHOL[$sp] } = split /\s/, $first_line; # ??? A MODIF
	# print "apres first_line ".__LINE__."\n";
	# close LNUP_LNWHOL;
	# ???
    

    if($bool_fixed_motif_threshold)
    {
	$rapports_display[$sp] = $rapports[$sp] = $motif_threshold;
    }
    else
    {
	$rapports[$sp] = $multipl_factor_motif_tresh * $upstr_nt[$sp] / 2 / $total_nt[$sp];
	$rapports_display[$sp] = sprintf "%2.2f ",$rapports[$sp];
    }

    ($not_verbose)or($sp==0)or print "rapports @rapports, rapports_display @rapports_display\n";
    print "$prog_tag A motif is interesting when its number of occurrences in merged\n";
    print "$prog_tag upsteam sequences is $multipl_factor_motif_tresh higher than in\n";
    print "$prog_tag the rest of the genome\n";
    print "$prog_tag Upstream seq cumulative length for $sp:$upstr_nt[$sp] (/ 2 because both orientations)\n";
    print "$prog_tag Total    seq cumulative length for $sp:$total_nt[$sp]\n";
    print "$prog_tag COMPUTED RATIO                 for $sp:$rapports_display[$sp]\n";    

    # creation of index files for annotation files (for faster searching and getting) ********************
    $index_annot_f[$sp] = $DATA_DIR.'INDEX_ANNOT_'."$name_fic_sp[$sp]".".txt";
    if(not -e $index_annot_f[$sp]){
	&create_index_file($way_to_annot_dir."$name_fic_sp[$sp]".".txt", $index_annot_f[$sp]);
    }
    #            (name of index files, ref to table storing byte positions of every annotations)
    &load_index_f($index_annot_f[$sp], \@{$annot_byte_position[$sp]} );
    print "index $index_annot_f[$sp] created\n";
    # print "VERIF annot_byte_position\n";
    # foreach(@{$annot_byte_position[$sp]}){ print "$_\n";}
    # print "END VERIF annot_byte_position\n";
    # exit;

    {
	no strict "refs";
	open(FANNOT.$sp, "$way_to_annot_dir"."$name_fic_sp[$sp]".".txt")or die "$prog_tag [Error] Cannot open $way_to_annot_dir"."$name_fic_sp[$sp]".".txt file line" .__LINE__."\n";
	
    }
} # end of for
if(not $not_verbose){
  print "verif handle_ref_glob\n";
  print "handle_ref_glob de 0 0: $handle_ref_glob00\n";
  print "handle_ref_glob de 0 1: $handle_ref_glob01\n";
  print "handle_ref_glob de 1 0: $handle_ref_glob10\n";
  print "handle_ref_glob de 1 1: $handle_ref_glob11\n";
  print "END verif handle_ref_glob\n";
  # exit;
  # open(FANNOT0, $way_to_annot_dir."$name_fic_sp[0]".".txt")or die "Can not open ".$way_to_annot_dir."$name_fic_sp[0]".".txt file line ".__LINE__."\n";
  # open(FANNOT1, $way_to_annot_dir."$name_fic_sp[1]".".txt")or die "Can not open ".$way_to_annot_dir."$name_fic_sp[1]".".txt file line ".__LINE__."\n";
  
  print "ID and SEQ files recorded, nb_whole_seq @nb_whole_seq, nb_upstr_seq @nb_upstr_seq, name_fic_sp @name_fic_sp\n";
}



# GLOBAL VARIABLES for MAIN prog
my %spacer_word1_word2_fic       = (); # hach table with five keys which allow to associate two related words
                                       # and file where the pair is found (keys

my @important_motifs             = (); # table of motifs by order of length and importance
my %H_important_motifs           = (); # hash which is the continuity of the table

# tie(%H_important_motifs, 'DB_File', 'H_important_motifs.txt', O_RDWR|O_CREAT,0666,$DB_BTREE);

# indice of the current size of @important_motifs
my $ind_important_motifs = 0;


my %trinuc1_trinuc2_spacer       = (); # hach table which associate two words of 3 nucleotides (one in each box) 
                                       # and a spacer. It contains a table of pointers on an input of the previous 
                                       # hach table
# tie %trinuc1_trinuc2_spacer, 'Tie::RefHash'; 
# tie %trinuc1_trinuc2_spacer, 'Tie::RefHash::Nestable'; 


my %fic_trinuc1_trinuc2          = (); # H table which associates a file to all its trinuc dyads
my %pos_mot                      = (); # table de hachage stockant les positions du mot
                                       # pour une sequence donnee
my %ID_diff_pos                  = (); # table de hachage qui enregistre les identifiant
                                       # des sequences concernees par un meme espacement entre deux boites 
my %words_for_a_given_spacer     = (); # table de hachage qui enregistre un comptage pour chaque couple (mots,espace
                                       # entre les positions d'un meme mot dans nos deux sequence)
my %exist_motif_diff_pos         = (); # table de hachage qui va permettre d'eviter les redondances de sous-mots
my %exist_spacer_word1_word2_fic = (); # table de hachage qui va permettre d'eviter les redondances de sous mots 
                                       # dans la table de hachage generale commune a tous nos fichiers    
my %nb_motifs_distincts_par_diff_par_fic = (); # cette table de hachage stocke le nombre d'occurrences (dans des 
                                       # sequences dostinctes  pour chaque difference de position possible dans les 
                                       # sequences comparees


# GLOBAL VAR used in RES_COMMON_ORTHO and called subprograms in
my $sp                   = 0;  # var giving species being treated
my $fic_sp;                    # var storing output file for good species
my $length_same_begin_init;    # const used to store length for getting subseq (see just before &RECUP...FIRST_CALL)
my @deb_prev_intit = ('zzzzzzzzzzzzzzzzzzzzzzzz')x2; # to know if we have to write again treatment of nnn, nnn...   

my %prev_infos_div       = (); # remember data used by get2.. subprog to avoid redondance in seq of each sp
my $entete               = ''; # string before results when we write in species files
my $intitules_ind_prem_l = ''; # record for the concerned line (number= key) the title

my @sorted_ind           = ();
my @seq_sp               = (); # table of parts of sp sequences
my @ID_f                 = (); # table of file IDs
my @infos_div            = (); # record data needed by get2.. subprog to avoid redondance in seq of each sp
my @pos_tri              = (); # record pos of trinuc -35 sp1 -10 sp1 -35 sp2 -10 sp2
my @pos_tri_tmp          = (); # save pos of letter for trinuc extensions
my @tab_lettre           = (); # pour recuperation des lettres (pour l'extension)

my $LOCAL_WIDTH_OF_L_HBOXES;   # length of trinuc of -35 box
my $LOCAL_WIDTH_OF_R_HBOXES;   # length of trinuc of -10 box

my @display              = (); # table used to record matrices associated with motifs
my %H_nbw                = (); # H table used to group and sort trinuc dyads by number of occurrences (in fact, 
                               # same H table as trinuc1_trinuc2_spacer (one key more, not really redundont H 
                               # table 
my @initial_border_words = (); # store initial border words for extension (useful for regexp with jokers for 
                               # MATRICE subprogramme)
my $alphabet_trinuc                 = 'acgt'; # used to count nb of fixed letter using a regexp on initial_border_words
my $initial_total_nb_boxes_letters  = 0;      # store nb of fixed letters compute on initial_border_words to compute during extension number of fixed letters
my $initial_total_boxes_length      = 0;      # store length of initial boxes before extension (used to compute number of fixed letters during extension
my $nb_of_l_after_extend_to_evaluate = 0;     # store the number of fixed letters in our two trinuc used for grouping		
			      
my @bool_first_call;           # initial boolean for extension when we call recup_lettre_sort (intialized just before)
                               # (used also in sort)


my $max_valid_prob_ext               = 0.15;  # max proba which had to be used for extension (print of warning if one 
                                             # of the proba of extension used is HIGHER to generate an interesting motif)
my @tempo_spacer_bounds_init = (100,0);      # permet d'initialiser les limites pour chaque lettre a la position 
                                             # d'interet pour calculer la longueur du spacer de l'expression reguliere
my @spacer_bounds        = ();
my @interv_int           = (); # will record limits of the intervals of sequences which are interesting according to 
                               # the middle of sequences scores. Only these intervals will be displayed
my $nb_mini_of_seq_by_interv = 4; # important var: fix minimal number of sequences used to extend motif (and so to compute proba)


# GLOBAL VAR used in VERIF_WORD_SPACE
my @words_to_sort_for_1pos       = ();

# GLOBAL VAR used in RECHERCHE
my %pos_mot_pour_ce_mot          = ();    
my %diff_trouvees_pour_ce_motif  = (); # utilise dans CALCUL_DIFF, 


# *********************************************************
# TREATMENT OF numeric directories (corresponding to biological functions)
# IN OTHER CAS, global treatment, we let go on

my $bool_dir = 0; # if 1 means it has got directories of directories

opendir(DIRS,$directory_for_ortho_seq) or die "$prog_tag [Error] Cannot open $directory_for_ortho_seq dir line ".__LINE__."\n";
my @dirs = readdir DIRS;
closedir DIRS;



# die "@dirs";
# my $STORE_DIR_strict = $directory_for_ortho_seq.'STORE_DIR'; # directory to store stacks for multiple low proba extension in sub SORT
# if(-d $STORE_DIR_strict){
    # print "on fait: rm -f $STORE_DIR_strict*\n";
   # system("rm -Rf $STORE_DIR_strict");
# }
# else{
  #  print "we create $STORE_DIR_strict\n";
# }
# $STORE_DIR_strict .= '/';

# mkdir($STORE_DIR_strict);

# my $STORE_DIR = $STORE_DIR_strict.'store_sub_SORT_stack';

# verif we have only directories and No_existing... file
DIR:for my $dir(@dirs)
{   # my $flag = 1;
    chomp($dir);

    #print "dir vaut $dir, ";
    #print "est-ce un repertoire?:";
    #(-d $dir)or next;
    #print "oui\n";
    #($dir =~ /^\.+/)and next;
    # if((! -d $dir)or($dir =~/^\.+/)){ 
    #print "dir $dir est un repertoire qui ne commence pas par un point\n";
    #   next DIR;	
    # }
    my $chem_dir = $directory_for_ortho_seq.$dir;
    print "\ndir : $dir; chem_dir= $chem_dir\n";
    
    # il suffit qu'il y ait un fichier de sequences pour que l'on considere qu'il n'y a pas des sous-repertoires a traiter
    if ($dir =~/^SIGffRid_ortholog_files\d+$/){ 
      $bool_dir = 1;
      ($not_verbose)or print "bool_dir mis a $bool_dir pour $dir\n";
      print " \n bool is: $bool_dir pour $dir\n";
      last DIR;
      
    }
   
#    else
#     {
# 	(-d $chem_dir)or print "dir $chem_dir non rep\n";
# 	(-e $chem_dir)or print "n' existe pas $chem_dir\n";
# 	print "ne ressemble pas a fic normal de rep ".($dir !~ /^(?:No_existing_promot_seq|ls|align\w+|CROSS\w+)\.txt$/ )."\n";
# 	print "ne ressemble pas a rep normal ".($dir !~ /^(?:\d+|\.+)$/ )."\n";
# 	print "fait partie des cas de repertoire\n";
#     }
}

my $dir = '';
if($bool_dir)
{
  # ($not_verbose)or print "bool_dir $bool_dir\n";
  for(my $i_dir = 0; $i_dir <= $#dirs; $i_dir++)
    {
      $dir = $dirs[$i_dir];
      if($dir !~ /^SIGffRid_ortholog_files\d+$/){ next; }
	# if((! -d $dir)or($dir =~/^\.+/)){ next; }
       # ($dir =~ /No_existing_promot_seq\.txt|STORE_DIR|^\.+/)and next;
	# print "dir $dir deux part\n";
	&MAIN_LOCAL;
    }

}
else
{
  # die "PB: booldir $bool_dir, dir $dir\n";
  &MAIN_LOCAL;
}

for my $i(0..$#name_fic_sp)
{
    my $filean = 'FANNOT'.$i; 
    close($filean);
}

close($handle_ref_glob00) or die "$prog_tag [Error] Impossible to close FSORSEQ$_, qui correspond a $handle_ref_glob00:$!, line ".__LINE__."\n";
close($handle_ref_glob01) or die "$prog_tag [Error] Impossible to close FSORID$_, qui correspond a $handle_ref_glob01:$!, line ".__LINE__."\n";
close($handle_ref_glob10) or die "$prog_tag [Error] Impossible to close FSORSEQ$_, qui correspond a $handle_ref_glob10:$!, line ".__LINE__."\n";
close($handle_ref_glob11) or die "$prog_tag [Error] Impossible to close FSORID$_, qui correspond a $handle_ref_glob11:$!, line ".__LINE__."\n";
# print `free`;

@t_mots_surrep = ();
%mots_surrep = ();

# system("echo \"\" | mail -s \"prog finished\" $author_mail");
# exit; #this is the end of the program hopefully Haribol!!!!!!!
# exit 1;

### ------------------------------------------------------
### ***********************START OF ALL SUBS******************

# ***********************************RECORD (UPDATE) THE STRING OF BOOLEAN FOR COUNTING TRINUC OCCURRENCES**********
# string of boolean is recorded in $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ]
# $ind_new_occ                        : indice of the corresponding composite word in array, 
# GLOBAL VAR %trinuc1_trinuc2_spacer  : H table used to group couples of interesting boxes, 
# $num_subw1, 
# $num_subw2, 
# $r_pack_box_spacer_trinuc, 
# $r_pack_ind_important_motifs_moins1 : keys to access to trinuc1_trinuc2_spacer H table, 
sub RECORD_TRINUC($$$\$\$)
{
    my ( $ind_new_occ, $num_subw1, $num_subw2, $r_pack_box_spacer_trinuc, $r_pack_ind_important_motifs_moins1) = @_;

    my $bool_exists = 1;

    # print "RECORD_TRINUC traitement de $num_subw1 $num_subw2 ".unpack('U',$$r_pack_box_spacer_trinuc).' '.unpack('N',$$r_pack_ind_important_motifs_moins1)."\n";
    # ($not_verbose)or print "SUB RECORD_TRINUC @_\n";
    (exists $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ])or $bool_exists = 0;

    if( ( $bool_exists )and( length( $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ] ) > $ind_new_occ ) )
    {
	substr( $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ], $ind_new_occ, 1 )++;
	# ($not_verbose)or print "incrementation de la i $ind_new_occ lettre de $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ]\n";
    }
    else
    {
	if( $bool_exists )
	{
	    # ($not_verbose)or print "longueur de chaine ".$trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ]." vaut ".length($trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ])."\n";
	    $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ]  .= ('0')x( $ind_new_occ - length( $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ] ) ).'1';
	    # ($not_verbose)or print "ajout en fin de chaine a la pos $ind_new_occ de '1' : $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ]\n";
	}
	else
	{
	    $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ] = ('0')x( $ind_new_occ ).'1';
	    # ($not_verbose)or print "creation de chaine a la pos $ind_new_occ on met '1' : $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $$r_pack_box_spacer_trinuc }{ $$r_pack_ind_important_motifs_moins1 }[ 1 ]\n";
	}
    }
}
# *********************************************************************************************************



# sub used for specific sorts
sub numerique{ $a <=> $b }
sub unpack_numerique_U{ unpack('U',$a) <=> unpack('U',$b) }
sub unpack_numerique_N{ unpack('N',$a) <=> unpack('N',$b) }

# ************************************BEGIN SUB VERIF********************************************************
# SUB which records potentially interesting words (owning to a couple of conserved boxes) in global structures so as to remove others in local structures (used for local search)
# param 1: ref to %words_for_a_given_spacer
# param 2: first key for %words_for_a_given_spacer H table (car on ne va parcourir qu'une valeur de la premiere cle)
# param 3: variable giving minimal space between two boxes
# param 4: variable giving maximal space between two boxes
# param 5: ref to %pos_mot (H table which give access to anonym array) store position for a given word
# param 6: ref to the array enumerating every possible shifts in the spacers
# param 7: ref to exist_motif_diff_pos H table, which allow to know if a word at a given position has already been recorded
sub VERIF_WORDS_SPACE($$$\%$\$)
{
    my ($cle1_sub,$min,$max,$pack_id_fic_seq,$r_ind_important_motifs) = @_;
    
    # ($not_verbose)or print "SUB VERIF_WORDS_SPACE @_\n";
    my %pos_mots_compares           = (); # hachage associant une position a des mot (hachage de tableau anonyme) (stocke aussi la position dans l'autre seq)
    my @tab_pos_mots_compares       = (); # idem mais tableau de positions ordonnees dans la deuxieme seq
    my @tab_bool_pos_mots_compares  = (); # tableau de booleen pour savoir quel seront les positions a effacer
    my @indice_pos_mots_compares    = (); # stocke l indice de la position dans le tableau pointe par le hachage %pos_mots_compares 
        
    # list of computed differences according to the real difference and the allowed spacer shift
    my @diff_shift = (); 
    
    # for a given difference (according to the allowed variation in the spacer
    for my $variat(@spacer_shift)
    {
	push @diff_shift , ($cle1_sub + $variat);
    
	# $diff_shift[$#diff_shift] = difference
	# cle2 = mot
	# pos = position pour cette difference et ce mot
	# on cree le hachage (de tableau) faisant correspondre a une position, le mot associe 
	
	for my $cle2(sort longueur keys %{ $words_for_a_given_spacer{$diff_shift[$#diff_shift]} })
	{
	    my $length_cle2 = length($cle2);

	    for my $ipos(0..$#{$pos_mot{$diff_shift[$#diff_shift]}{$cle2}})
	    {
		my $pos = ${$pos_mot{$diff_shift[$#diff_shift]}{$cle2}}[ $ipos ];
		my $bool_not_inserted = 1;
		for my $i(0..$#{$pos_mots_compares{$pos}[0]})
		{
		    if( length( ${$pos_mots_compares{$pos}[0]}[$i] )< $length_cle2 )
		    { 
			$bool_not_inserted = 0; 
			@{$pos_mots_compares{$pos}[0]}[ $i+1..$#{$pos_mots_compares{$pos}[0]}+1 ] = @{$pos_mots_compares{$pos}[0]}[$i..$#{$pos_mots_compares{$pos}[0]}];
			@{$pos_mots_compares{$pos}[1]}[ $i+1..$#{$pos_mots_compares{$pos}[1]}+1 ] = @{$pos_mots_compares{$pos}[1]}[$i..$#{$pos_mots_compares{$pos}[1]}];
			@{$pos_mots_compares{$pos}[2]}[ $i+1..$#{$pos_mots_compares{$pos}[2]}+1 ] = @{$pos_mots_compares{$pos}[2]}[$i..$#{$pos_mots_compares{$pos}[2]}];
			${$pos_mots_compares{$pos}[0]}[ $i ] = $cle2;
			# variation of the spacer in the other seq
			${$pos_mots_compares{$pos}[1]}[ $i ] = $variat;
			# position of the motif in the other seq
			${$pos_mots_compares{$pos}[2]}[ $i ] = $pos - $diff_shift[$#diff_shift];
		    }
		}
		if( $bool_not_inserted )
		{ 
		    push @{$pos_mots_compares{$pos}[0]}, $cle2;
		    push @{$pos_mots_compares{$pos}[1]}, $variat;
		    push @{$pos_mots_compares{$pos}[2]}, $pos - $diff_shift[$#diff_shift];
		}
		# push @{$pos_mots_compares{$pos}[0]}, $cle2;
		# @{$pos_mots_compares{$pos}[0]} = reverse sort longueur @{$pos_mots_compares{$pos}[0]};
	    }
	}
    }
    
    # affichage de verification de la table de hachage %pos_pos_compares et de son contenu (tableau de mots)
    # for my $val(sort numerique keys(%pos_mots_compares)){
	# for my $chaque_pos(0..$#{$pos_mots_compares{$val}[0]}){
	 #   print "cle $val val ${$pos_mots_compares{$val}[0]}[$chaque_pos]\nvariat ${$pos_mots_compares{$val}[1]}[$chaque_pos] posmot2 ${$pos_mots_compares{$val}[2]}[$chaque_pos]\n"; 
	# }
    # }
    # print "fin H $cle1_sub\n\n";
    # die;
    # on cree le tableau des positions triees et celui des mots (dans le meme ordre)
    
    my $cptline_tab_bool_pos = 0;
    
    for my $pos_ordo( sort numerique keys %pos_mots_compares)
    {
	for my $chaque_pos(0..$#{$pos_mots_compares{$pos_ordo}[0]})
	{
	    push @tab_pos_mots_compares, $pos_ordo;
	    push @indice_pos_mots_compares, $chaque_pos;
	    $tab_bool_pos_mots_compares[$cptline_tab_bool_pos][0] = 0;
	    $tab_bool_pos_mots_compares[$cptline_tab_bool_pos++][1] = 0;
	}
    }
    
    # if(! $not_verbose)
    # {
	# print "affichage de controle\n";
	# for my $pos_ordo(0..$#tab_pos_mots_compares){
	   # print "pos $tab_pos_mots_compares[$pos_ordo] indice $indice_pos_mots_compares[$pos_ordo] mot ${$pos_mots_compares{$tab_pos_mots_compares[$pos_ordo]}[0]}[$indice_pos_mots_compares[$pos_ordo]]\n";
	# }
	# print "fin controle\n";
    # }

    for my $pos_ordo(0..$#tab_pos_mots_compares-1)
    {
	my $first_word = ${ $pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[0] }[ $indice_pos_mots_compares[$pos_ordo] ];
	my $length_first_word = length($first_word);
	
	POS2:for my $pos_ordo2($pos_ordo+1..$#tab_pos_mots_compares)
	{
	    # ON NE PEUT PAS GENERALISER UN CALCUL EN SE BASANT UNIQUEMENT SUR LES POSITIONS, CAR LA DISTANCE ENTRE LES DEUX BOITES DEPEND AUSSI DE LA LONGUEUR DU MOT (VARIABLE), MAIS ON DOIT POUVOIR ACCELERER LES OPERATIONS
	    
	    my $box_spacer = $tab_pos_mots_compares[$pos_ordo2]-$tab_pos_mots_compares[$pos_ordo]-$length_first_word;
	    
	    last POS2 if($box_spacer > $max);
	    # si la position qui suit est inferieure a la position courante + ecart mini entre les boites, rien ne servira de faire la comparaison, on fait donc avancer le compteur de debut de boucle (accelere le traitement)
	    next POS2 if($box_spacer < $min);
	    
	    my $pack_box_spacer = pack('U',$box_spacer);
	    $box_spacer = ();

	    my $second_word = ${$pos_mots_compares{$tab_pos_mots_compares[$pos_ordo2]}[0]}[$indice_pos_mots_compares[$pos_ordo2]];
	    my $length_second_word = length($second_word);

	    if( (!exists $exist_spacer_word1_word2_fic{ $pack_box_spacer }{ $first_word }{ $second_word }{ pack('n',$tab_pos_mots_compares[$pos_ordo]) } )and
		((! $bool_use_filter_for_mini_nb_of_comm_letters)or
		 ($length_first_word + $length_second_word >= $mini_number_of_letters_in2words)) )
	    {
		if(($length_second_word + $length_first_word < $mini_number_of_letters_in2words)and($bool_use_filter_for_mini_nb_of_comm_letters)){ die "$prog_tag [Error] PB filtre, line ".__LINE__."\n"; }
		# delete $exist_spacer_word1_word2_fic{ $pack_box_spacer }{ $first_word }{ $second_word }{ pack('n',$tab_pos_mots_compares[$pos_ordo])};
		
		# the first case of tab_bool_pos... is a boolean we put to 1 if the motif corresponds to a -35 box
		if(! $tab_bool_pos_mots_compares[$pos_ordo][0])
		{  
		    $tab_bool_pos_mots_compares[$pos_ordo][0] = 1; 
		}
		
		# the second case of tab_bool_pos... is a boolean we put to 1 if the motif corresponds to a -10 box
		if(! $tab_bool_pos_mots_compares[$pos_ordo2][1])
		{  
		    $tab_bool_pos_mots_compares[$pos_ordo2][1] = 1; 
		}
		
		# record of interesting motifs (no redondance) in a global table
		
		my $packed_ind_dyad = ();
		
		if( scalar(keys( %{$spacer_word1_word2_fic{ $pack_box_spacer }{ $first_word }{ $second_word }})) > 0)
		{
		    my ($k) = keys(%{$spacer_word1_word2_fic{ $pack_box_spacer }{ $first_word }{ $second_word }});
		    $packed_ind_dyad = $spacer_word1_word2_fic{ $pack_box_spacer }{ $first_word }{ $second_word }{ $k };

		    # if(($first_word =~ /gaa/)and($second_word =~ /gtt/))
		    # {
			# print "SIGR $first_word $second_word ".unpack('U',$pack_box_spacer)." FILE ".$tab_ID_fic[unpack('n', $pack_id_fic_seq)][0]."\n";
		    # }
		    # ($not_verbose)or print "k (idficseq) trouve vaut ".unpack('n',$k)." correspondant a ".unpack('N',$packed_ind_dyad)."\n";
		    # ($not_verbose)or print "we find an already recorded word in file ".$tab_ID_fic[unpack('n', $pack_id_fic_seq)][0]."\n first_word $first_word second_word $second_word pack_box_spacer ".unpack('U',$pack_box_spacer)." must correspond to ".$r_important_motifs->[unpack('N',$packed_ind_dyad)][0].' '.unpack('U',$r_important_motifs->[unpack('N',$packed_ind_dyad)][1]).' '.$r_important_motifs->[unpack('N',$packed_ind_dyad)][2]." GET EXISTING\n";
		    # record of file id

#ICI
# if($useDiskStorage)
# {

# }
# else
# {
		    push @{$H_important_motifs{ $packed_ind_dyad }[0] }, $pack_id_fic_seq;
		    # record of the position of -35 box in this file for the first bacteria
		    push @{$H_important_motifs{ $packed_ind_dyad }[1] }, pack('n',$tab_pos_mots_compares[$pos_ordo]);

                    # ($not_verbose)or print "we record packed pos1 ".$tab_pos_mots_compares[$pos_ordo]." GET EXISTING\n";

		    # record of the position of -35 box in this file for the other bacteria
		    push @{$H_important_motifs{ $packed_ind_dyad }[2] }, pack('n',${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[2]}[ $indice_pos_mots_compares[$pos_ordo] ]);
		    # ($not_verbose)or print "we record packed pos2 ".${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[2]}[ $indice_pos_mots_compares[$pos_ordo] ]." GET EXISTING\n";
		    # spacer shift in the 2 seq
		    push @{$H_important_motifs{ $packed_ind_dyad }[3] }, pack('cc',${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[1]}[ $indice_pos_mots_compares[$pos_ordo] ],${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo2] }[1]}[ $indice_pos_mots_compares[$pos_ordo2] ]);

		    # ($not_verbose)or print "we record packed shifts ".${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[1]}[ $indice_pos_mots_compares[$pos_ordo] ].' '.${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo2] }[1]}[ $indice_pos_mots_compares[$pos_ordo2] ]." GET EXISTING\n";

		    # ($not_verbose)or print "on met $tab_pos_mots_compares[$pos_ordo] dans r_H_important_motifs ind $packed_ind_dyad a la pos $#{$H_important_motifs{ $packed_ind_dyad }[1] }\n";
		    # record in a hash of usefull characteristics not to record twice the same motif
		    # first word
		    # second word
		    # box spacer
		    # ID of the file
		    $spacer_word1_word2_fic{ $pack_box_spacer }{ $first_word }{ $second_word }{ $pack_id_fic_seq } = $packed_ind_dyad;
# }		    
		    # to avoid results for the same position and motif
		}
		else
		{
		    my $pack_ind_important_motifs = pack('N',$$r_ind_important_motifs);
		    $important_motifs[$$r_ind_important_motifs][0] = $first_word;
		    $important_motifs[$$r_ind_important_motifs][1] = $pack_box_spacer;
		    $important_motifs[$$r_ind_important_motifs][2] = $second_word;
		   
		    # ($not_verbose)or print "add in important_motifs table $first_word, ".unpack('U',$pack_box_spacer).", $second_word for file ".$tab_ID_fic[unpack('n',$pack_id_fic_seq)][0]." shifts ".${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[1]}[ $indice_pos_mots_compares[$pos_ordo] ].' '.${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo2] }[1]}[ $indice_pos_mots_compares[$pos_ordo2] ]." positions $tab_pos_mots_compares[$pos_ordo] ".${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[2]}[ $indice_pos_mots_compares[$pos_ordo] ]." file ".$tab_ID_fic[unpack('n',$pack_id_fic_seq)][0].' line '.__LINE__."\n"; 
		    push @{$H_important_motifs{$pack_ind_important_motifs}[0]}, $pack_id_fic_seq;
		    push @{$H_important_motifs{$pack_ind_important_motifs}[1]}, pack('n',$tab_pos_mots_compares[$pos_ordo]);
		    push @{$H_important_motifs{$pack_ind_important_motifs}[2]}, pack('n',${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[2]}[ $indice_pos_mots_compares[$pos_ordo] ]);
		    push @{$H_important_motifs{$pack_ind_important_motifs}[3]}, pack('cc',${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo] }[1]}[ $indice_pos_mots_compares[$pos_ordo] ],${$pos_mots_compares{ $tab_pos_mots_compares[$pos_ordo2] }[1]}[ $indice_pos_mots_compares[$pos_ordo2] ]);
		    
		    # record in a hash of usefull characteristics not to record twice the same motif
		    # first word
		    # second word
		    # box spacer
		    # ID of the file)
		    $spacer_word1_word2_fic{ $pack_box_spacer }{ $first_word }{ $second_word }{ $pack_id_fic_seq } = $pack_ind_important_motifs;
		    
		    $$r_ind_important_motifs++;
		}
		
		my $pack_ind_important_motifs_moins1;
		if(! defined $packed_ind_dyad )
		{          
		    my $a_empaqueter = $$r_ind_important_motifs - 1;
		    $pack_ind_important_motifs_moins1 = pack('N',$a_empaqueter);
		}
		else
		{
		    $pack_ind_important_motifs_moins1 = $packed_ind_dyad;
		}
		
		# if(! $not_verbose)
		# {
		    # my $unpack_i = unpack('N',$pack_ind_important_motifs_moins1)-1;
		    # print "prev ind_important_motifs corresponds to motif ".$important_motifs[$unpack_i][0].' '.unpack('U',$important_motifs[$unpack_i][1]).' '.$important_motifs[$unpack_i][0]."\n";
		# }

		# ENREGISTREMENT POUR LE RECOUPEMENT GENERAL DES RESULTATS
		# pour les positions allant de 0 a la longueur du mot moins la longueur du sous-mot regarde (3)
		
		# JUSQUE LA RIEN NE CHANGE POUR TRINUC A GAP ???

		# GROS CHANGEMENTS *********************************************************************

		
		# $trinuc_lengths->[0] = minimal length of trinuc
		# for each length of trinuc (sorted by length)
		for my $id_trinuc_l1(0..$#{$trinuc_lengths[0]})
		{
		    # every_keys_for_first_word_and_this_trinuc_length
		    for my $ind_basic_re( $trinuc_lengths[1][$id_trinuc_l1][0]..$trinuc_lengths[1][$id_trinuc_l1][1] )
		    {
			# de 0 a la longueur du mot moins la taille du premier trinuc de regroupement
			for my $diff_pos_sousmots(0..($length_first_word - $trinuc_lengths[0][$id_trinuc_l1]))
			{
			   # print "on cherche $tmp_VAR_trinuc[$ind_basic_re] de longueur $trinuc_lengths[0][$id_trinuc_l1] dans\n".substr($first_word,$diff_pos_sousmots,$trinuc_lengths[0][$id_trinuc_l1])."\nextrait de $first_word de la pos $diff_pos_sousmots de longueur $trinuc_lengths[0][$id_trinuc_l1]\net on veut faire un remplacement dans $VAR_trinuc[$ind_basic_re], ind_basic $ind_basic_re correspond a $tmp_VAR_trinuc[$ind_basic_re]\n";

			    substr($first_word,$diff_pos_sousmots,$trinuc_lengths[0][$id_trinuc_l1]) =~ /$tmp_VAR_trinuc[$ind_basic_re]/;
			    my $first_key = $VAR_trinuc[$ind_basic_re];

			    # print "first_key $first_key devient \n";

			    my ($dol1,$dol2,$dol3);
			    if(defined $1){ $dol1 = $1; }else{ $dol1 = (); die "$prog_tag [Error] dol1 not defined line ".__LINE__."\n"; }
			    if(defined $2){ $dol2 = $2; }else{ $dol2 = (); }
			    if(defined $3){ $dol3 = $3; }else{ $dol3 = (); }
			    if(defined $dol1){ $first_key =~ s/w+/$dol1/; }
			    if(defined $dol2){ $first_key =~ s/w+/$dol2/; }
			    if(defined $dol3){ $first_key =~ s/w+/$dol3/; }

			    # print "first_key $first_key\n";

			    # if(length($first_key) != $trinuc_lengths[0][$id_trinuc_l1]){ die "line ".__LINE__."\n"; }
			    # if(length($first_key) != length($VAR_trinuc[$ind_basic_re])){ die "line ".__LINE__."\n"; }
			    
			    for my $id_trinuc_l2(0..$#{$trinuc_lengths[0]})
			    {
				SECTRI:for my $ind_basic_re2( $trinuc_lengths[1][$id_trinuc_l2][0]..$trinuc_lengths[1][$id_trinuc_l2][1] )
				{
				    # if(($bool_use_filter_for_mini_nb_of_comm_letters)and($mini_number_of_letters_in2trinuc > $nbl[$ind_basic_re] + $nbl[$ind_basic_re2] )){
					# print "prochain car bool_use_filter_for_mini_nb_of_comm_letters $bool_use_filter_for_mini_nb_of_comm_letters and mini_number_of_letters_in2trinuc $mini_number_of_letters_in2trinuc > nbl de $ind_basic_re $nbl[$ind_basic_re] + nbl de $ind_basic_re2 $nbl[$ind_basic_re2]\n";
				     # }
				     # else{
					# print "PAS prochain car bool_use_filter_for_mini_nb_of_comm_letters $bool_use_filter_for_mini_nb_of_comm_letters and mini_number_of_letters_in2trinuc $mini_number_of_letters_in2trinuc > nbl de $ind_basic_re $nbl[$ind_basic_re] + nbl de $ind_basic_re2 $nbl[$ind_basic_re2]\nUne condition non respectee\n";
					# exit;
				     # }

				    # ($bool_use_filter_for_mini_nb_of_comm_letters)and($mini_number_of_letters_in2trinuc > ($nbl[$ind_basic_re] + $nbl[$ind_basic_re2]) )and next;
				    if($bool_use_filter_for_mini_nb_of_comm_letters)
				    {
					my $cptl2t = $nbl[$ind_basic_re] + $nbl[$ind_basic_re2];
					(($mini_number_of_letters_in2trinuc > $cptl2t )or($maxi_number_of_letters_in2trinuc < $cptl2t ))and next;
				    }

				    # de 0 a la taille du deuxieme trinuc de regroupement
				    for my $diff_pos_sousmots2(0..($length_second_word - $trinuc_lengths[0][$id_trinuc_l2]))
				    {
					# print "on cherche $tmp_VAR_trinuc[$ind_basic_re2] de longueur $trinuc_lengths[0][$id_trinuc_l2] dans\n".substr($second_word,$diff_pos_sousmots2,$trinuc_lengths[0][$id_trinuc_l2])."\nextrait de $second_word de la pos $diff_pos_sousmots2 de longueur $trinuc_lengths[0][$id_trinuc_l2]\net on veut faire un remplacement dans $VAR_trinuc[$ind_basic_re2], ind_basic $ind_basic_re2 correspond a $tmp_VAR_trinuc[$ind_basic_re2]\n";
					my $sec_key = $VAR_trinuc[$ind_basic_re2];	
					$not_verbose or print " chk 1 $first_key and $sec_key \n";
					if($first_key =~ /^[^n]{3}$/ and $sec_key eq 'www' and $bool_www_treated){
							next SECTRI;
					}
					
					substr($second_word,$diff_pos_sousmots2,$trinuc_lengths[0][$id_trinuc_l2]) =~ /$tmp_VAR_trinuc[$ind_basic_re2]/;
					# print "sec_key here $sec_key devient\n";
					if(defined $1){ $dol1 = $1; }else{ $dol1 = (); }
					if(defined $2){ $dol2 = $2; }else{ $dol2 = (); }
					if(defined $3){ $dol3 = $3; }else{ $dol3 = (); }
					if(defined $dol1){ $sec_key =~ s/w+/$dol1/; }
					if(defined $dol2){ $sec_key =~ s/w+/$dol2/; }
					if(defined $dol3){ $sec_key =~ s/w+/$dol3/; }

					# print "sec_key $sec_key\n";
					# print "sec_key $sec_key car tmp_VAR_trinuc de $ind_basic_re2 $tmp_VAR_trinuc[$ind_basic_re2]\nsubstr vaut ".substr($second_word,$diff_pos_sousmots2,$trinuc_lengths[0][$id_trinuc_l2])."\nVAR_trinuc de $ind_basic_re2 $VAR_trinuc[$ind_basic_re2]\n";

					 if(length($sec_key) != length($VAR_trinuc[$ind_basic_re2])){ die "$prog_tag [Error] difference of length!! line ".__LINE__."\n"; }
					# die "tmp_VAR_trinuc de $ind_basic_re2 $VAR_trinuc[$ind_basic_re2]\n";

					# EN COURS
					
					# my $subw1 = substr($first_word,$diff_pos_sousmots,3);
					# my $subw2 = substr($second_word,$diff_pos_sousmots2,3);
					
					# my $num_subw1 = $H_trinuc{ substr($first_word,$diff_pos_sousmots,$trinuc_l) };
					my $num_subw1 = $H_trinuc{ $first_key };
					my $num_subw2 = $H_trinuc{ $sec_key };
					# print "first_key $first_key correspond a $H_trinuc{ $first_key } dans la Htable donc a $table_trinuc[$H_trinuc{ $first_key }]\n";
					
					# my $num_subw2 = $H_trinuc{ substr($second_word,$diff_pos_sousmots2,$trinuc_l) };
					# pack_box_spacer_trinuc = length between words 
					# + position of second trinuc in second word (diff_pos_sousmots2)
					# + length of first word
					# - length of first trinuc ($trinuc_lengths[0][$id_trinuc_l1])
					# - position of first trinuc in first word (diff_pos_sousmots)
					my $pack_box_spacer_trinuc = pack('U',unpack('U',$pack_box_spacer)+$diff_pos_sousmots2+$length_first_word-$diff_pos_sousmots-$trinuc_lengths[0][$id_trinuc_l1]);
					
					my $bool_not_already_found = -1;
					
					# print "NOUVEAU $table_trinuc[$num_subw1] $table_trinuc[$num_subw2] $pack_box_spacer_trinuc\n";
					# print "test existence num_subw1 $num_subw1 ";
					# print "val $table_trinuc[$num_subw1] ";
					# print "num_subw2 $num_subw2 ";
					# print "val $table_trinuc[$num_subw2] ";
					# print "pack_box_spacer_trinuc ".(unpack('U',$pack_box_spacer)+$diff_pos_sousmots2+$length_first_word-$diff_pos_sousmots-$trinuc_lengths[0][$id_trinuc_l1])."\n";

					if((! defined $num_subw2)or(! defined $num_subw1)){ exit; }
					if(exists $trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }{ $pack_ind_important_motifs_moins1 }[ 0 ])
					{
					    for my $idt(0..$#{$trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }{ $pack_ind_important_motifs_moins1 }[ 0 ] })
					    {
						if(unpack('U',${$trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }{ $pack_ind_important_motifs_moins1 }[ 0 ] }[$idt]) ==  $diff_pos_sousmots)
						{
						    $bool_not_already_found = $idt;
						    last;
						}
					    }
					}
					
					if($bool_not_already_found == -1)
					{
					    push @{$trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }{ $pack_ind_important_motifs_moins1 }[ 0 ]},pack('U',$diff_pos_sousmots);
					    # ${$trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $box_spacer_trinuc }{ $pack_ind_important_motifs_moins1 }[ 1 ]}[ $#{$trinuc1_trinuc2_spacer{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }{ $pack_ind_important_motifs_moins1 }[ 0 ]} ][ $#{$H_important_motifs{$pack_ind_important_motifs_moins1}[1]} ] = 1;
					    # ($not_verbose)or print "RECORD_TRINUC for $table_trinuc[$num_subw1] $table_trinuc[$num_subw2] ".unpack('U',$pack_box_spacer_trinuc)." file ".$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_ind_important_motifs_moins1}[0]}[$#{$H_important_motifs{$pack_ind_important_motifs_moins1}[1]}])][0]."\n";
					    &RECORD_TRINUC( $#{$H_important_motifs{$pack_ind_important_motifs_moins1}[1]}, $num_subw1, $num_subw2, \$pack_box_spacer_trinuc, \$pack_ind_important_motifs_moins1 );
					    
					    # hash table to fusion words overlapped in a same sequence (first value = num of implicated word in global table, second value = indice of the involved position in the table pointed by the global hash table ( $important_motifs{}
#			    push @{$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ pack('U',$box_spacer_trinuc) }[0]}, pack('N',$pack_ind_important_motifs_moins1);		
					    # print "on ajoute pour le fichier ".$tab_ID_fic[unpack('n',$pack_id_fic_seq)][0]." pour les sous-mots $table_trinuc[$num_subw1] ".unpack('U',$pack_box_spacer_trinuc)." $table_trinuc[$num_subw2] ln 604 le mot ".$important_motifs[unpack('N',$pack_ind_important_motifs_moins1)][0].' '.unpack('U',$important_motifs[unpack('N',$pack_ind_important_motifs_moins1)][1]).' '.$important_motifs[unpack('N',$pack_ind_important_motifs_moins1)][2]."\n";	
					    push @{$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }[0]}, $pack_ind_important_motifs_moins1;			
					    push @{$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }[1]}, pack('n',$#{$H_important_motifs{$pack_ind_important_motifs_moins1}[1]});
					    
					}
					else
					{
					    my $bool_not_fic_record = -1;
					    for my $idf(reverse 0..$#{$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }[0]} )
					    {
						if( ( unpack('N',${$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }[0]}[ $idf ]) == unpack('N',$pack_ind_important_motifs_moins1) )and( unpack('n',${$H_important_motifs{$pack_ind_important_motifs_moins1}[1]}[ unpack('n',${$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }[1]}[ $idf ]) ]) == $tab_pos_mots_compares[$pos_ordo] ) )
						{
						    $bool_not_fic_record = unpack('n',${$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }[1]}[ $idf ]);
						    last;
						}
					    }
					    
					    
					    
					    if( $bool_not_fic_record == -1 )
					    {
						# print "on ajoute pour le fichier ".$tab_ID_fic[unpack('n',$pack_id_fic_seq)][0]." pour les sous-mots $table_trinuc[$num_subw1] ".unpack('U',$pack_box_spacer_trinuc)." $table_trinuc[$num_subw2] ln 625 le mot ".$important_motifs[unpack('N',$pack_ind_important_motifs_moins1)][0].' '.unpack('U',$important_motifs[unpack('N',$pack_ind_important_motifs_moins1)][1]).' '.$important_motifs[unpack('N',$pack_ind_important_motifs_moins1)][2]."\n";	
						push @{$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }[0]}, $pack_ind_important_motifs_moins1;
						push @{$fic_trinuc1_trinuc2{ $pack_id_fic_seq }{ $num_subw1 }{ $num_subw2 }{ $pack_box_spacer_trinuc }[1]}, pack('n',$#{$H_important_motifs{$pack_ind_important_motifs_moins1}[1]});
						
				# ($not_verbose)or print "RECORD_TRINUC for $table_trinuc[$num_subw1] $table_trinuc[$num_subw2] ".unpack('U',$pack_box_spacer_trinuc)." file ".$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_ind_important_motifs_moins1}[0]}[$#{$H_important_motifs{$pack_ind_important_motifs_moins1}[1]}])][0]."\n";
						&RECORD_TRINUC( $#{$H_important_motifs{$pack_ind_important_motifs_moins1}[1]}, $num_subw1, $num_subw2, \$pack_box_spacer_trinuc, \$pack_ind_important_motifs_moins1 );
					    }
					    else
					    {
						# ($not_verbose)or print "RECORD_TRINUC for $table_trinuc[$num_subw1] $table_trinuc[$num_subw2] ".unpack('U',$pack_box_spacer_trinuc)." file ".$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_ind_important_motifs_moins1}[0]}[ $bool_not_fic_record ])][0]."\n";
						&RECORD_TRINUC( $bool_not_fic_record, $num_subw1, $num_subw2, \$pack_box_spacer_trinuc, \$pack_ind_important_motifs_moins1 );
					    }
					    
					}
				    }
				} # AJOUT
			    } # AJOUT
			} # AJOUT			
		    }
		} # AJOUT
		
		# FIN?? GROS CHANGEMENTS *********************************************************************

		# hachage enorme (mais efface a chaque nouveau fichier) permettant d'eviter les redondances de deux boites dues a des sous mots dans notre hachage general visant a croiser les resultats entre fichiers)
		for my $longueur(($ordre_mini_moins1 + 2) .. $length_first_word)
		{
		    for my $longueur2(($ordre_mini_moins1 + 2) .. $length_second_word)
		    {
			# pour les positions allant de 0 a la longueur du mot moins la longueur du sous-mot regarde
			for my $diff_pos_sousmots(0..($length_first_word - $longueur))
			{
			    for my $diff_pos_sousmots2(0..($length_second_word - $longueur2))
			    {
				# print "on eteint mot1 ".substr($first_word,$diff_pos_sousmots,$longueur).",boxspacer ".(unpack('U',$pack_box_spacer)+$diff_pos_sousmots2+$length_first_word-$diff_pos_sousmots-$longueur).", mot2 ".substr( $second_word,$diff_pos_sousmots2,$longueur2 )."\n";
				# keys: spacer, subword1, subword2, pos
				$exist_spacer_word1_word2_fic{ pack('U',unpack('U',$pack_box_spacer)+$diff_pos_sousmots2+$length_first_word-$diff_pos_sousmots-$longueur) }{ substr($first_word,$diff_pos_sousmots,$longueur) }{ substr( $second_word,$diff_pos_sousmots2,$longueur2)}{ pack('n',$tab_pos_mots_compares[$pos_ordo]+$diff_pos_sousmots) } = 1;
			    }
			}
		    }
		}
	    }
	}
    }
    for my $pos_a_eff(reverse 0..$#tab_bool_pos_mots_compares)
    {
	if((!$tab_bool_pos_mots_compares[$pos_a_eff][0])and(!$tab_bool_pos_mots_compares[$pos_a_eff][1]))
	{
	    delete $words_for_a_given_spacer{$cle1_sub}{${$pos_mots_compares{$tab_pos_mots_compares[$pos_a_eff]}[0]}[$indice_pos_mots_compares[$pos_a_eff]]}; # ??
	    $#{$pos_mot{$cle1_sub}{${$pos_mots_compares{$tab_pos_mots_compares[$pos_a_eff]}[0]}[$indice_pos_mots_compares[$pos_a_eff]] }} = -1;
	    delete $pos_mot{$cle1_sub}{${$pos_mots_compares{$tab_pos_mots_compares[$pos_a_eff]}[0]}[$indice_pos_mots_compares[$pos_a_eff]]};
	    delete $pos_mot{$cle1_sub}{${$pos_mots_compares{$tab_pos_mots_compares[$pos_a_eff]}[1]}[$indice_pos_mots_compares[$pos_a_eff]]};
	    delete ${$pos_mots_compares{$tab_pos_mots_compares[$pos_a_eff]}[0]}[$indice_pos_mots_compares[$pos_a_eff]];
	    delete ${$pos_mots_compares{$tab_pos_mots_compares[$pos_a_eff]}[1]}[$indice_pos_mots_compares[$pos_a_eff]];
	    delete ${$pos_mots_compares{$tab_pos_mots_compares[$pos_a_eff]}[2]}[$indice_pos_mots_compares[$pos_a_eff]];
	    
	    if( $#{ $pos_mots_compares{ $tab_pos_mots_compares[$pos_a_eff]}[0]} == -1)
	    { 
		delete $pos_mots_compares{ $tab_pos_mots_compares[$pos_a_eff] }; 
	    }
	
	    splice(@tab_bool_pos_mots_compares,$pos_a_eff,1);
	    splice(@tab_pos_mots_compares,$pos_a_eff,1);
	}
    }
    
    # ICI TRI EN FONCTION DE LA LONGUEUR (traitement des motifs longueur apres longueur (sous prog pour pouvoir eviter recherche des sous_mots
    # ******************************************?? VERIFIER QU ON A DEJA TOUTES LES DIFFERENCES POUR UN MOT DONNE
    
    # partie permettant (par appel) de preciser que le tri se fait sur la longueur des chaines
    # sub longueur{ length($a) <=> length($b) }
    
    
    @words_to_sort_for_1pos = ();
    my $prev_pos = 'a';
    
    # on commence par trier en fonction dela position car on ne supprimera un mot que s il a la meme position
    # qu un mot plus grand faisant partie de la meme boite

    # for all the positions of the words sorted in increasing order
    for my $pos(sort numerique keys(%pos_mots_compares))
    {
	
	# for all the word corresponding to this position
	for my $ind_tab_anonyme(0..$#{$pos_mots_compares{$pos}[0]})
	{
	    (!defined ${$pos_mots_compares{$pos}[0]}[$ind_tab_anonyme])or push @words_to_sort_for_1pos, ${$pos_mots_compares{$pos}[0]}[$ind_tab_anonyme];
	}
	
	# treatment once we have got all the words for a given position (so, treatment is shifted)
	# sort of the words according to their length
	@words_to_sort_for_1pos = reverse sort longueur @words_to_sort_for_1pos;
	
	# suppression des occurrences de mots eventuellement indesirables
	&TREATMENT_REMOVE_SUBWORDS($pos,$cle1_sub);
	
	# apres traitement, nous vidons le tableau
	$#words_to_sort_for_1pos = -1;
	$prev_pos = $pos;
    }
    
    @words_to_sort_for_1pos = reverse sort longueur @words_to_sort_for_1pos;
    
    &TREATMENT_REMOVE_SUBWORDS($prev_pos,$cle1_sub);
    @words_to_sort_for_1pos = ();
}
# ************************************END SUB VERIF**********************************************************

# ************************************MAIN *****************************************************************

sub MAIN_LOCAL
{
    # INITIALIZATIONS

    # ($not_verbose)or print "SUB MAIN @_\n"; 
    # hach table with five keys which allow to associate two related words and file where the pair is found (keys 
    %spacer_word1_word2_fic = ();
    # table of motifs by order of length and importance
    @important_motifs = ();
    %H_important_motifs = (); # hash which is the continuity of the table
    # indice of the current size of @important_motifs
    $ind_important_motifs = 0;

    # hach table which associate two words of 3 nucleotides (one in each box) and a spacer. It contains a table of pointers on an input of the previous hach table
    %trinuc1_trinuc2_spacer = ();
    # tie %trinuc1_trinuc2_spacer, 'Tie::RefHash';
    # tie %trinuc1_trinuc2_spacer, 'Tie::RefHash::Nestable'; 

    %fic_trinuc1_trinuc2 = ();

    
    %pos_mot       = ();        # table de hachage stockant les positions du mot pour une sequence donnee
    %ID_diff_pos   = ();        # table de hachage qui enregistre les identifiants des sequences concernees par un meme espacement entre deux boites 
    %words_for_a_given_spacer     = (); # table de hachage qui enregistre un comptage pour chaque couple (mots,espace entre les positions d'un meme mot dans nos deux sequence)
    %exist_motif_diff_pos         = (); # table de hachage qui va permettre d'eviter les redondances de sous-mots
    %exist_spacer_word1_word2_fic = (); # table de hachage qui va permettre d'eviter les redondances de sous mots dans la table de hachage generale commune a tous nos fichiers    

    chomp $directory_for_ortho_seq;
    chomp $dir;
    my $current_chemin = Cwd::cwd();
#     print "\n\n my current_chemin is: $current_chemin \n\n";
    # 	my $ch=<STDIN>;
    
    
    if ($bool_dir) {chdir("$directory_for_ortho_seq$dir")or die "$prog_tag [Error] chdir $directory_for_ortho_seq$dir impossible: $!, line ".__LINE__."\n";
			print "\n\n my current_chemin is: $directory_for_ortho_seq$dir \n\n";
   			  #	my $ch=<STDIN>;
    
    }
    else {chdir("$directory_for_ortho_seq")or die "$prog_tag [Error] chdir $directory_for_ortho_seq impossible: $!, line ".__LINE__."\n";
    }

    # ($not_verbose)or print "treatment of $directory_for_ortho_seq$dir directory\n";
    print "Treatment of directory: $directory_for_ortho_seq$dir \n";

    # open(PWDDebug, "pwd |")or die "Impossible to open ls command line ".__LINE__.": $!";
    # print "Nous sommes dans ".<PWDDebug>."\n";
    # close PWDDebug;

    # ?? DEMANDER POUR SUPPRIMER SI UN ls.txt EXISTE
    # recherche des noms des fichiers d'orthologues dans le repertoire sense les contenir
   # remove qw(ls.txt* CROSS*.txt No_existing_promot_seq.txt*);
    unlink glob("ls.txt*");
    unlink glob("CROSS*.txt");
    unlink glob("No_existing_promot_seq.txt*");
    unlink($align_id1, $align_id2);

   # system('rm -f ls.txt*');
   # system("rm -f CROSS*".'.txt*');
   # system("rm -f $align_id1");
   # system("rm -f $align_id2");
   # system("rm -f No_existing_promot_seq.txt*");

    # system("rm -f ls.txt CROSS$name_fic_sp[0]"."_$name_fic_sp[1]".'.txt'." align$name_fic_sp[0]".'.txt'." align$name_fic_sp[1]".'.txt');

    # system('ls > ls.txt')or die "Impossible to do ls line ".__LINE__.": $!";
   # open(FIC_ORTHO,'ls |') or die "Impossible to do ls: $!";
    my @fic_ortho = ();
    opendir(FIC_ORTHO, '.')or die "$prog_tag [Error] Impossible to open current directory line ".__LINE__.": $!\n";
    @fic_ortho = grep -T,  readdir FIC_ORTHO;
    closedir FIC_ORTHO;
    
    if($bool_pc_actif)
    {
	# open(NB,'wc -l ls.txt |')or die "Impossible to count number of lines for ls.txt file, line ".__LINE__.": $!";

	# determination de la valeur seuil (nombre de fichiers concernes par un motif pour considerer qu'un motif est interessant)
	# my $nbl = <NB>;
	# $nbl =~ /^\s*(\d+)/;
	$min_nb_of_waited_sigma_motifs = $min_nb_of_waited_sigma_motifs * $#fic_ortho;
	# close NB or die "Impossible to close NB: $!";
    }


    # record length of upstream (or promot seq) used foreach file and bacteria to compute real position in results
    my @length_upstr_seq = ();

    # pour chacun de ces fichiers de sequences amonts d'orthologues
    for my $fic_seq_all_promot(@fic_ortho)
    {
	chomp $fic_seq_all_promot;

# 	my $ch=<STDIN>;
	if( ($fic_seq_all_promot =~ /^(?:align|ls\.txt)/) or ($fic_seq_all_promot =~ /^res_sort.*\.txt/) 
	      or ($fic_seq_all_promot =~/^result_grep.*\.txt/)or (-d $fic_seq_all_promot) ){ next; }
	if ($fic_seq_all_promot =~ /No_existing_promot_seq\.txt|STORE_DIR/){ next};
	print " Now treating: $fic_seq_all_promot\n";
	# (-d $fic_seq_all_promot)and next;

	# ($not_verbose)or print "we treat $fic_seq_all_promot file\n";

	# we record file ID for results display
	# push @tab_ID_fic, $fic_seq_all_promot;

	%pos_mot       = (); # table de hachage stockant les positions du mot pour une sequence donnee
	%ID_diff_pos   = (); # table de hachage qui enregistre les identifiants des sequences concernees par un meme espacement entre deux boites
	%words_for_a_given_spacer = (); # table de hachage qui enregistre un comptage pour chaque couple (mots,espace entre les positions d'un meme mot dans nos deux sequence)
	%exist_motif_diff_pos = ();  # table de hachage qui va permettre d'eviter les redondances de sous-mots
	%exist_spacer_word1_word2_fic = (); # table de hachage qui va permettre d'eviter les redondances de sous mots dans la table de hachage generale commune a tous nos fichiers

	# on l'ouvre
	open(FPROM,"< $fic_seq_all_promot")or die "$prog_tag [Error] Impossible to open $fic_seq_all_promot: $!, line ".__LINE__."\n";

	# creation du fichier de sortie dans le meme repertoire que les fichiers de sequences d'orthologues
	(!$bool_edit_fic_paire_ortho) or open(FSOR,"> $entete_fic_res_1_pair".$fic_seq_all_promot) or die "$prog_tag [Error] impossible de creer un fichier qui s'appelle $entete_fic_res_1_pair".$fic_seq_all_promot.": $!, line ".__LINE__."\n";

	%nb_motifs_distincts_par_diff_par_fic = (); # cette table de hachage stocke le nombre d'occurrences (dans des sequences dostinctes  pour chaque difference de position possible dans les sequences comparees

	# va permettre de connaitre le rang du motif pour un ordre donne
	my $classement_du_motif_par_ordre = 0;
	my $length_prec_word = $ordre_maxi + 2;

	# store length of seq
	my @l_seq = ();
	$#t_mots_surrep == -1 and die "$prog_tag [Error] t_mots_surrep empty!!!: @t_mots_surrep, line ".__LINE__."\n";
	for my $reg_exp_MM(@t_mots_surrep)
	{
	    if($length_prec_word == length($reg_exp_MM))
	    {
		$classement_du_motif_par_ordre++;
	    }
	    else{
		$classement_du_motif_par_ordre = 1;
	    }

	    &RECHERCHE($reg_exp_MM, $classement_du_motif_par_ordre, \@l_seq);

	}
	# we record file ID for results display, and length of each seq
	# print "$fic_seq_all_promot, @l_seq\n";
	if($#l_seq == -1){ die "$prog_tag [Error] longueur de seq non determinee dans $fic_seq_all_promot!!, line ".__LINE__."\n"; }
	push @tab_ID_fic, [$fic_seq_all_promot, @l_seq];
	# push @length_upstr_seq, [@l_seq];

	# cle1 = difference
	# pour toutes les differences...
	for my $cle1(sort numerique keys %nb_motifs_distincts_par_diff_par_fic)
	{
	    my $verif = 0;

	    # on compte le nombre de motifs distincts repertories en considerant la variabilite possible de l espacement entre boites
	    for my $varia(@spacer_shift)
	    {
		my $cle1_plus_shift = $cle1 + $varia;
		if(exists $nb_motifs_distincts_par_diff_par_fic{$cle1_plus_shift})
		{
		    $verif += $nb_motifs_distincts_par_diff_par_fic{$cle1_plus_shift};
		}
	    }

	    # si le nombre de motifs repertories pour cette difference est superieur a 1, cela va nous interesser
	    if($verif >= $nb_seq_min_owning_a_word_in_one_ortholog_file)
	    {
		&VERIF_WORDS_SPACE($cle1,$diff_min,$diff_max,pack('n',$#tab_ID_fic),\$ind_important_motifs);

		# cle2 = mot
		for my $cle2(reverse sort longueur keys %{ $words_for_a_given_spacer{$cle1} })
		{
		    my $nb_occu_word_for_intervalle = 0;
		    for my $varia(@spacer_shift)
		    {
			my $cle1_plus_shift = $cle1 + $varia;
			if(exists $words_for_a_given_spacer{$cle1_plus_shift}{$cle2})
			{
			    $nb_occu_word_for_intervalle += $words_for_a_given_spacer{$cle1_plus_shift}{$cle2};
			}
		    }

		    if($nb_occu_word_for_intervalle > 1)
		    {
			if($bool_edit_fic_paire_ortho)
			{
			    for my $varia(@spacer_shift)
			    {
				my $cle1_plus_shift = $cle1 + $varia;
				if(exists $words_for_a_given_spacer{$cle1_plus_shift}{$cle2})
				{
				    # print "nombre de mots sup a 1 $words_for_a_given_spacer{$cle1}{$cle2}\n";
				    print FSOR "diff $cle1_plus_shift mot $cle2 compte $words_for_a_given_spacer{$cle1_plus_shift}{$cle2} pos @{$pos_mot{$cle1_plus_shift}{$cle2}} IDseq $ID_diff_pos{$cle1_plus_shift}{$cle2}\n";
				    # ($not_verbose)or print "diff $cle1_plus_shift mot $cle2 compte $words_for_a_given_spacer{$cle1_plus_shift}{$cle2} pos @{$pos_mot{$cle1_plus_shift}{$cle2}[0]} IDseq $ID_diff_pos{$cle1_plus_shift}{$cle2}\n";
				}
			    }
			}
		    }
		    else
		    {
			# ($not_verbose)or print "(1)efface diff $cle1 mot $cle2 compte $words_for_a_given_spacer{$cle1}{$cle2}\n";
			delete $ID_diff_pos{$cle1}{$cle2};
			delete $words_for_a_given_spacer{$cle1}{$cle2};
			$#{$pos_mot{$cle1}{$cle2}} = -1;
			delete $pos_mot{$cle1}{$cle2};
		    }
		}
		# &CAPTURE_COMPUTER_STAT("MEM LOAD file\n");
	    }
	    # elsif(!$not_verbose)
	    # {
		# print "le nombre de motifs concernant la difference $cle1 est inferieur a 1 ($verif)\n";
	    # }

	}
	close FPROM;

	if($bool_edit_fic_paire_ortho)
	{
	    close FSOR;
	    print "Le fichier genere s'appelle $1".'motifs_'."$fic_seq_all_promot\n";
	}
    }

    # &CAPTURE_COMPUTER_STAT("MEM just before EXTEND_DYAD\n");

    print "Beginning of words fusion for a same file (EXTEND_DYAD)\n";

    # le placer ailleurs (dans autre boucle pour diminuer complexite et memoire vive necessaire)
    &EXTEND_DYAD(\$ind_important_motifs);

    %fic_trinuc1_trinuc2    = ();
    %spacer_word1_word2_fic = ();

    # print "size main table (important_motifs) ".($#important_motifs+1)."\n";

    # &CAPTURE_COMPUTER_STAT("MEM JUST BEFORE CROSSING RESULTS\n");

    print "End of words fusions for a same file (EXTEND_DYAD), fic_trinuc1_trinuc2 and spacer_word1_word2 H removed\n";
    print "Compute of number of files involved in a trinucleotids dyad set, and display of results\n";


    # creating ID file and seq file where searching interesting regexp_motifs from fasta files corresponding to upstream seq

    # ICI
    &RES_COMMON_ORTHO;

    #untie %trinuc1_trinuc2_spacer;
    %trinuc1_trinuc2_spacer = ();
    @important_motifs       = ();
    %H_important_motifs     = ();

    print "results file written\n";


   # close FIC_ORTHO;
   # system("rm -f ls.txt*");
    unlink glob("ls.txt*");

    my $fin_sigffrid = time();
    my $infos_last = '';

	my $duree = $fin_sigffrid - $debut_sigffrid;
    if($duree >= 0)
    {
	$infos_last = "The programme began at: ".(scalar localtime($debut_sigffrid))."\nThe programme stopped at: ".(scalar localtime($fin_sigffrid))."\nThe programme last $duree secunds\n";
	print $infos_last;
    }
    else{
		print join("\n", 	"Problem computer clock was modified while program was running",
							"duree:$infos_last",
							"debut:$debut_sigffrid",
							"fin:$fin_sigffrid",
							"line ".__LINE__."\n");
	}

    if($bool_print_cross)
    {
	print FCROSS $infos_last;
	close(FCROSS);
    }

    chdir("$current_chemin")or die "$prog_tag [Error] chdir $current_chemin impossible: $!, line ".__LINE__."\n";

}
# ************************************END MAIN *************************************************************

# OK gap
# SUB REMOVE SUBWORDS*********************************************************************
sub TREATMENT_REMOVE_SUBWORDS($$\$)
{
    # reference to the list of all the words at a given position
    # position of the current word
    # reference to reference to a hach table used not to record subwords contained by a good larger box occurrence
    # difference which is the same for all the treated word in this subprogram (if we are taking account of the variable length of the spacer)
    # reference to reference to a table of all the allowed spacer shifts
    #
    #
    my ($pos_mot_en_cours_sub, $difference) = @_;

    # ($not_verbose)or print "SUB TREATMENT_REMOVE_SUBWORDS @_\n";

    # cle2 = treated word
    for my $cle2(@words_to_sort_for_1pos)
    {
	
	my $bool_repertorie = 0; # boolean which indicates if a word is already recorded in another relashionship
	
	# for each one of the possible spacer shifts
      REPER:for my $spacer_shift(@spacer_shift)
      {
	  
	  # we compute the real difference according to the allowed shift in the spacer
	  my $difference_shift = $difference + $spacer_shift;
	  
	  # if the word is already recorded with the sames parameters, we affect 1 to our boolean
	  if(exists $exist_motif_diff_pos{$cle2}{$difference_shift}{$pos_mot_en_cours_sub})
	  {
	      # ($$ID_diff_pos{$difference_shift}{$cle2} =~ /$gene_IDs_sub/)){ # et dans la meme sequence
	      $bool_repertorie = 1;
	      # and we go to the end of the loop as we don't need to verify following different spacer shifts
	      last REPER; 
	  }
      }
	
	# si ce mot n'a pas ete trouve dans un autre motif a la meme position 
	if($bool_repertorie)
	{
	    # suppression de l occurrence de ce mot puisqu elle est deja presente a l interieur d un autre mot a la meme position pour une difference similaire
	    # print "on supprime $cle2 car deja present\n";
	    delete $ID_diff_pos{$difference}{$cle2};
	    delete $words_for_a_given_spacer{$difference}{$cle2};
	    $#{$pos_mot{$difference}{$cle2}} = -1;
	    delete $pos_mot{$difference}{$cle2};
	}
	else
	{
	    my $length_cle2 = length($cle2);
	    
	    # list of ccomputed differences according to the real difference and the allowed spacer shift
	    my @diff_shift = (); 
	    
	    for my $diff_var(@spacer_shift)
	    {
		push @diff_shift , ($difference + $diff_var);
	    }
	    
	    # pour eviter redondance de sous-motifs
	    # pour les longueurs allant de la plus petite longueur admise par l'ordre du modele de markov jusqu'a la longueur de notre mot moins 1
	    for my $longueur(($ordre_mini_moins1 + 2) .. ($length_cle2))
	    {
		# pour les positions allant de 0 a la longueur du mot moins la longueur du sous-mot regarde
		for my $diff_pos_sousmots(0..($length_cle2 - $longueur))
		{
		    for my $diff_spacer(@diff_shift)
		    {
			$exist_motif_diff_pos{substr($cle2,$diff_pos_sousmots,$longueur)}{$diff_spacer}{$pos_mot_en_cours_sub + $diff_pos_sousmots} = 1;
		    }
		}
	    }
	}
    }	  	  
}
# END SUB REMOBE SUBWORDS****************************************************************



# sub programme used in EXTEND_DYAD to remove old hash table records and create new ones 
sub REMOVE($$$$$$\$\$\$\$\@\@\@\@\$\$)
{
    my($local_pack_pos_reelle_deb, 
       $local_pack_pos_reelle_deb_autre_seq, 
       $local_pack_shiftsp_autre_seq, 
       $local_ind_deb,
       $local_i,
       $local_diff_w1,
       $num_w1,
       $num_w2,
       $r_pack_file,
       $r_pack_boxsp,
       $r_champs_num,
       $r_champs_ind_diff_pos,
       $r_champs_ind_fic_tri, 
       $r_champs_pos_reelle, 
       $r_chaine1,
       $r_ind_important_motifs) = @_;
    
    # ($not_verbose)or print "SUB REMOVE @_\n";

#    if(($table_trinuc[$num_w1] eq 'cga')and($table_trinuc[$num_w2] eq 'cgc')and(unpack('U',$$r_pack_boxsp) == 21))
#    {
#	print "\n\n***********A L'ORIGINE***********\n";
#	for my $i_supp(reverse $local_ind_deb..$local_i)
#	{
#	    my $id_big_table = $r_champs_num->[$i_supp];
	    
#	    if(exists $$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $id_big_table }[ 1 ])	
#	    {
#		
#		print "on a rr_trinuc1_trinuc2_spacer $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp)." $id_big_table 1 de $r_champs_ind_diff_pos->[$i_supp] indHmotif ".unpack('n',${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ])." qui est a ".substr($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $id_big_table }[ 1 ],unpack('n',${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ]) ,1).": $important_motifs[ $r_champs_num->[$i_supp] ][ 0 ] $important_motifs[ $r_champs_num->[$i_supp] ][ 1 ] $important_motifs[ $r_champs_num->[$i_supp] ][ 2 ]\n";
#	    }
#	}
#	print "***************************\n";
#    }
    # record of new word
    $$r_chaine1 =~ /(.+?)([-]+)(.+?)$/;
    my $mot1 = undef;

    if(defined $1)
    {
	$mot1 = $1;
    }
    
    if((not defined $1)or(not defined $2))
    {
	print join("\n", "$prog_tag [REMOVE] parameters:",
		   "local_pack_pos_reelle_deb:".unpack('n',$local_pack_pos_reelle_deb), 
		   "local_pack_pos_reelle_deb_autre_seq:".unpack('n', $local_pack_pos_reelle_deb_autre_seq), 
		   "local_pack_shiftsp_autre_seq:".unpack('cc',$local_pack_shiftsp_autre_seq), 
		   "local_ind_deb:$local_ind_deb",
		   "local_i:$local_i",
		   "local_diff_w1:$local_diff_w1",
		   "num_w1:$num_w1 corresponding to $table_trinuc[$num_w1] seed",
		   "num_w2:$num_w2 corresponding to $table_trinuc[$num_w2] seed",
		   "r_pack_file:".unpack('n',$$r_pack_file),
		   "r_pack_boxsp:".unpack('U',$$r_pack_boxsp),
		   "r_champs_num:".unpack('N',$r_champs_num->[$local_ind_deb]),
		   "r_champs_ind_diff_pos:".$r_champs_ind_diff_pos->[$local_ind_deb],
		   "r_champs_ind_fic_tri:".join(',', @$r_champs_ind_fic_tri), 
		   "r_champs_pos_reelle:".join(',',@$r_champs_pos_reelle), 
		   "r_chaine1:".$$r_chaine1,
		   "r_ind_important_motifs:".$$r_ind_important_motifs,
		   "line ".__LINE__."\n");
    }
    if(not defined $1)
    {
	warn "$prog_tag [REMOVE][warn] dol1 not defined for r_chaine1:".$$r_chaine1." using regexp '(.+?)([-]+)(.+?)': seeds merged line ".__LINE__."\n";
    }
    if(not defined $2)
    {
	warn "$prog_tag [REMOVE][warn] dol2 not defined for r_chaine1:".$$r_chaine1." using regexp '(.+?)([-]+)(.+?)': seeds merged line ".__LINE__."\n";
    }
    my $pack_espace = undef; # pack('U',length($2));
    my $mot2 = undef;

    if(!defined $3)
    { 
	$pack_espace = pack('U',0);
	
	# $mot2 = '-'; # initially, changed 2019 02 25

	# arbitrary separation 2019 02 25
	# merged seeds
	# not good, ideally, we must know how to distinguish the two seeds
	# to set clearly mot1 and mot2 and deal with motifs under the form:
	# mot1_n{0,5}_mot2 with a kind of optional spacer
	# this initialisation if dol1 not defined was made in 2018 02 26

	my $verif_table_trinuc_num_w1 = $table_trinuc[$num_w1];
	my $verif_table_trinuc_num_w2 = $table_trinuc[$num_w2];
	&REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w1);
	&REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w2);
	# mot1 has the longest possible chain attributed
	if( $$r_chaine1 =~ /^(.*?$verif_table_trinuc_num_w1.*?)($verif_table_trinuc_num_w2.*?)$/ )
	{
	    $mot1 = $1;
	    $mot2 = $2;
	    print "$prog_tag we record mot1:$mot1, mot2:$mot2 for chain $$r_chaine1 that match second seed $table_trinuc[$num_w2], line ".__LINE__."\n";
	}
	else
	{
	    die "$prog_tag [Error] regexp ^(.*?$table_trinuc[$num_w1].*?)($table_trinuc[$num_w2].*?)\$ with first seed $table_trinuc[$num_w1] and second seed $table_trinuc[$num_w2] not found in chaine1 $$r_chaine1, line ".__LINE__."\n";
	}
    } 
    else{
	$pack_espace = pack('U',length($2));	
	$mot2 = $3; 
    }

    #  ($not_verbose)or STDERR print "chaine1 $$r_chaine1 dans remove\n";
    my $pack_indice_big_words = (); # record the indice of the big dyad word if it is already recorded so as not to remove it
  
    if(! $not_verbose)
    {
	my $verif_table_trinuc_num_w1 = $table_trinuc[$num_w1];
	my $verif_table_trinuc_num_w2 = $table_trinuc[$num_w2];
	&REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w1);
	&REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w2);
	if($mot1 !~ /$verif_table_trinuc_num_w1/){ die "$prog_tag [Error] PB: mot1 $mot1 does not contain rw1 $table_trinuc[$num_w1], line ".__LINE__."\n"; }
	if(($mot2 ne '-')and($mot2 !~ /$verif_table_trinuc_num_w2/)){ die "$prog_tag [Error] PB: mot2 $mot2 does not contain rw2 $table_trinuc[$num_w2], line ".__LINE__."\n"; }
	
	if($#{ $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $r_champs_num->[$local_ind_deb] }[ 0 ] }== -1){ die "$prog_tag [Error] RIEN file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." rw1 $table_trinuc[$num_w1] rw2 $table_trinuc[$num_w2] boxsp ".unpack('U',$$r_pack_boxsp)." mot1 $important_motifs[ unpack('N',$r_champs_num->[$local_ind_deb]) ][0] ".unpack('U',$important_motifs[ unpack('N',$r_champs_num->[$local_ind_deb]) ][1])." $important_motifs[ unpack('N',$r_champs_num->[$local_ind_deb]) ][2]\nlocalind deb $local_ind_deb locali $local_i, line ".__LINE__."\n";}
    }
    # ($not_verbose)or print "FILE ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]."\n";

    my $pack_diff_pos = ${ $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $r_champs_num->[$local_ind_deb] }[ 0 ] }[ $r_champs_ind_diff_pos->[$local_ind_deb] ];

    if(!defined $pack_diff_pos){ die "$prog_tag [Error] diff_pos non defini pour w1  $table_trinuc[$num_w1] w2 $table_trinuc[$num_w2] file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." sp ".unpack('U',$$r_pack_boxsp)." localinddeb ".unpack('N',$r_champs_num->[$local_ind_deb])." file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0].", line ".__LINE__."\n"; }   

    my $bool_not_existing_diff_pos = -1; # will record eventually the indice of the word if it already exists (with all the sames characteristics

    # die "boxsp ".unpack('U',$$r_pack_boxsp)."\n";

    # IF NEW WORD DOES NOT EXIST IN THE GLOBAL TABLE, WE CREATE IT OK
    if(scalar(keys(%{$spacer_word1_word2_fic{ $pack_espace }{ $mot1 }{ $mot2 }}))==0 )
    {
	my $pack_n_r_ind_important_motifs = pack('N',$$r_ind_important_motifs);
	$important_motifs[$$r_ind_important_motifs][0] = $mot1;
	$important_motifs[$$r_ind_important_motifs][1] = $pack_espace;
	$important_motifs[$$r_ind_important_motifs][2] = $mot2;
	
	push @{$H_important_motifs{$pack_n_r_ind_important_motifs}[0]}, $$r_pack_file;
	push @{$H_important_motifs{$pack_n_r_ind_important_motifs}[1]}, $local_pack_pos_reelle_deb;
	push @{$H_important_motifs{$pack_n_r_ind_important_motifs}[2]}, $local_pack_pos_reelle_deb_autre_seq;
	push @{$H_important_motifs{$pack_n_r_ind_important_motifs}[3]}, $local_pack_shiftsp_autre_seq;
	# ($not_verbose)or print "We record new word in REMOVE: $mot1 ".unpack('U', $pack_espace)." $mot2, indice $$r_ind_important_motifs positions ".unpack('n', $local_pack_pos_reelle_deb).' '.unpack('n', $local_pack_pos_reelle_deb_autre_seq)." shiftsp ".unpack('cc', $local_pack_shiftsp_autre_seq).", file ".$tab_ID_fic[unpack('n', $$r_pack_file)][0]."\n";

	$spacer_word1_word2_fic{ $pack_espace }{ $mot1 }{ $mot2 }{ $$r_pack_file } = $pack_n_r_ind_important_motifs;
	# print "on ajoute pour le fichier ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." pour les sous-mots $table_trinuc[$num_w1] ".unpack('U',$$r_pack_boxsp)." $table_trinuc[$num_w2] ln 1044 le mot ".$important_motifs[unpack('N',$pack_n_r_ind_important_motifs)][0].' '.unpack('U',$important_motifs[unpack('N',$pack_n_r_ind_important_motifs)][1]).' '.$important_motifs[unpack('N',$pack_n_r_ind_important_motifs)][2]."\n";	

	push @{$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 0 ]}, $pack_n_r_ind_important_motifs;
	push @{$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ]}, pack('n',$#{$H_important_motifs{ $pack_n_r_ind_important_motifs }[0]});
	
	push @{ $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_n_r_ind_important_motifs }[ 0 ] }, $pack_diff_pos;     
	$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_n_r_ind_important_motifs }[ 1 ] = (('0')x( $#{$H_important_motifs{ $pack_n_r_ind_important_motifs }[0]})).'1';     

	$pack_indice_big_words = $pack_n_r_ind_important_motifs;

	# ($not_verbose)or print "PACK_INDICE_BIG_WORDS ".unpack('N',$pack_indice_big_words)." NEW WORD\n\n";
	$$r_ind_important_motifs++;
    }
    else
    {
	my $bool_not_fic_diff = 1;

	# IF THE WORD EXISTS BUT NOT FOR THE SAME FILE

	if(!exists $spacer_word1_word2_fic{ $pack_espace }{ $mot1 }{ $mot2 }{ $$r_pack_file })
	{
	    # we get the indice of the related word already recorded 
	    for my $f(keys %{ $spacer_word1_word2_fic{ $pack_espace }{ $mot1 }{ $mot2 }})
	    {
		$pack_indice_big_words = $spacer_word1_word2_fic{ $pack_espace }{ $mot1 }{ $mot2 }{ $f };
		# ($not_verbose)or print "PACK_INDICE_BIG_WORDS ".unpack('N',$pack_indice_big_words)."WE RECORD AN ALREADY EXISTING WORD (BUT RELATED UNTIL THIS MOMENT ONLY TO ANOTHER OR OTHER FILE(S)), ITS ID IS ".unpack('N',$pack_indice_big_words)." file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]."\n";

		$bool_not_fic_diff = 0;
		last;
	    }	
	    if( $bool_not_fic_diff )
	    {
		if(!defined $pack_indice_big_words){  print "Indice_big_words not defined\n"; }
		die "$prog_tag [Error] PB: H for space ".unpack('U',$pack_espace)." word1 $mot1 word2 $mot2 no existing in ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." file but we cannot get any other ID related!!, line ".__LINE__."\n";
	    }
	}
	else
	{
	    $pack_indice_big_words = $spacer_word1_word2_fic{ $pack_espace }{ $mot1 }{ $mot2 }{ $$r_pack_file };
	    # ($not_verbose)or print "PACK_INDICE_BIG_WORDS ".unpack('N',$pack_indice_big_words)."WORD ALREADY RECORDED FOR THIS FILE ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]."\n";
	    if(!defined $pack_indice_big_words){  print "pack_indice_big_words not defined\n"; }
	}

	my $bool_not_already_recorded  = 1;

	# CE MOT EXISTE DEJA POUR CE MEME FICHIER
	if($bool_not_fic_diff)
	{
	    # if there is already a record related to this file, we search for eventually the same position
	    # ($not_verbose)or print "on met a jour mot1 $mot1 espace ".unpack('U',$pack_espace)." mot2 $mot2 pour le fichier ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]."\n";

	    # booleen pour se memoriser si l'on a trouve la position coorespondante ou non
	    my $bool_real_pos_found = -1;

	    for my $ind_champs_pos_reelle($local_ind_deb..$local_i)
	    {

		if( $r_champs_pos_reelle->[ $ind_champs_pos_reelle ] != unpack('n',$local_pack_pos_reelle_deb) )
		{
		    if( $bool_real_pos_found != -1)
		    {
			# la position reelle a deja ete trouvee, aucun enregistrement correspondant dans trinuc n'a pu etre trouve, nous en creons donc 1
			push @{$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[0]}, $pack_diff_pos;
			
			# ($not_verbose)or print "RECORD_TRINUC for $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp)." file ".$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_indice_big_words}[0]}[ unpack('n', ${$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[1] }[ $r_champs_ind_fic_tri->[ $bool_real_pos_found ] ]) ])][0]."\n";
			&RECORD_TRINUC( unpack('n', ${$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[1] }[ $r_champs_ind_fic_tri->[ $bool_real_pos_found ] ]), $num_w1, $num_w2, $r_pack_boxsp, \$pack_indice_big_words );

			# ( $not_verbose )or print "WE RECORD A NEW TRINUC, FOR AN EXISTING FICTRI file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." 1er\n";

			$bool_not_already_recorded = 0;
			if($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[1] !~ /^\d+$/){ die "$prog_tag [Error] line ".__LINE__." non numerique\n";}
			last;
		    }
		   # else
		   # {
			# la position reelle n'a pas ete trouvee pour le moment
			# ($not_verbose)or print "NO FIC_TRINUC1_TRINUC2 DEFINED FOR THIS WORD, SO NOT_ALREADY_RECORDED maybe file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]."\n";
		   # }
		}
		elsif( unpack('N',$r_champs_num->[ $ind_champs_pos_reelle ]) == unpack('N',$pack_indice_big_words) )
		{
		    # ($not_verbose)or print "FIC_TRINUC1_TRINUC2 DEFINED FOR THIS WORD\n";
		
		    $bool_real_pos_found = $ind_champs_pos_reelle;
			    
		    # we search for eventual same diff pos for the same position
		   
		    if( unpack('U',${$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[ 0 ]}[ $r_champs_ind_diff_pos->[ $ind_champs_pos_reelle ] ]) == unpack('U', $pack_diff_pos) )
		    {
			$bool_not_existing_diff_pos = $ind_champs_pos_reelle;
			# case where diff_pos is the same, sequence is already taken in account
			# ($not_verbose)or print "A RECORD WITH SAME ID, FILE AND POSITION WAS FOUND IN H_IMPORTANT_MOTIFS file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]."\n";

			# if we find the same position, we know that some data are already recorded
			$bool_not_already_recorded = 0;
			last;
		    }
		    elsif( $ind_champs_pos_reelle == $local_i )
		    {
			# la position reelle a deja ete trouvee, aucun enregistrement correspondant dans trinuc n'a pu etre trouve, nous en creons donc 1
			push @{$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[0]}, $pack_diff_pos;
			
			# ($not_verbose)or print "RECORD_TRINUC for $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp)." file ".$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_indice_big_words}[0]}[ unpack('n', ${$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[1] }[ $r_champs_ind_fic_tri->[ $ind_champs_pos_reelle ] ]) ])][0]."\n";
			&RECORD_TRINUC( unpack('n',${$fic_trinuc1_trinuc2{ $$r_pack_file }{ $H_trinuc{$num_w1} }{ $num_w2 }{ $$r_pack_boxsp }[1] }[ $r_champs_ind_fic_tri->[ $ind_champs_pos_reelle ] ]), $num_w1, $num_w2, $r_pack_boxsp, \$pack_indice_big_words );

			# ( $not_verbose )or print "WE RECORD A NEW TRINUC, FOR AN EXISTING FICTRI file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." 2eme\n";
			if($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[1] !~ /^\d+$/){ die "$prog_tag [Error] line ".__LINE__." non numerique\n";}
			$bool_not_already_recorded = 0;
			last;
		    }
		}
	    }
	}


	if($bool_not_already_recorded)
	{
	    
	    # ($not_verbose)or print "AUCUNE POSITION EQUIVALENTE N'A ETE TROUVEE DANS FIC_TRINUC1_TRINUC2\n";
	    # boolean to know if related trinuc record exists
	    my $bool_not_half_recorded = 1;
	    
	    if($#{$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[ 0 ]} == -1)
	    {
		# $bool_not_half_recorded = 1;
		# ($not_verbose)or print "NO RECORD FOR THIS WORD IN TRINUC1_TRINUC2_SPACER w1 $table_trinuc[$num_w1] w2 $table_trinuc[$num_w2] boxsp ".unpack('U',$$r_pack_boxsp)." idbigword ".unpack('N',$pack_indice_big_words)."\n";
	    }
	    else
	    {
		# if the same position exists, we get if it corresponds to the same diffpos for the motif (generally, 2 values max)
	      for my $loc_diff_pos(0..$#{$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[ 0 ]})
	      {
		  if( unpack('U',${$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[ 0 ]}[$loc_diff_pos]) == unpack('U',$pack_diff_pos) )
		  {
		      # ($not_verbose)or print "WORD RECORDED FOR THIS DIFFPOS IN TRINUC1_TRINUC2_SPACER w1 $table_trinuc[$num_w1] w2 $table_trinuc[$num_w2] boxsp ".unpack('U',$$r_pack_boxsp)." idbigword ".unpack('N',$pack_indice_big_words).", THIS WORD CORRESPONDS TO A DIFFERENT POSITION OF THE CURRENT\n";
		      # if diffpos is the same, we increase trinuc1_trinuc2_spacer only if it concerns a new file or if the position in the file for the already recorded result is different 
		      
		      push @{$H_important_motifs{$pack_indice_big_words}[0]}, $$r_pack_file;
		      push @{$H_important_motifs{$pack_indice_big_words}[1]}, $local_pack_pos_reelle_deb;
		      push @{$H_important_motifs{$pack_indice_big_words}[2]}, $local_pack_pos_reelle_deb_autre_seq;
		      push @{$H_important_motifs{$pack_indice_big_words}[3]}, $local_pack_shiftsp_autre_seq;
		      
		      # ($not_verbose)or print "word already recorded, file ".$tab_ID_fic[unpack('n', $$r_pack_file)][0].' '.$important_motifs[unpack('N',$pack_indice_big_words)][0].' '.unpack('U',$important_motifs[unpack('N',$pack_indice_big_words)][0]).' '.$important_motifs[unpack('N',$pack_indice_big_words)][2].', we add positions pos1 '.unpack('n',$local_pack_pos_reelle_deb).' pos2 '.unpack('n',$local_pack_pos_reelle_deb).' shifts '.unpack('cc',$local_pack_shiftsp_autre_seq)." line ".__LINE__."\n";
		      # ($not_verbose)or print "tableau important motifs, tailles: pour 0 $#{$H_important_motifs{$pack_indice_big_words}[0]} ,pour 1 $#{$H_important_motifs{$pack_indice_big_words}[1]} ,pour 2 $#{$H_important_motifs{$pack_indice_big_words}[2]} , pour 3 $#{$H_important_motifs{$pack_indice_big_words}[3]}\n";

		      # ($not_verbose)or print print "on ajoute pour le fichier ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." pour les sous-mots $table_trinuc[$num_w1] ".unpack('U',$$r_pack_boxsp)." $table_trinuc[$num_w2] ln 1189 le mot ".$important_motifs[unpack('N',$pack_indice_big_words)][0].' '.unpack('U',$important_motifs[unpack('N',$pack_indice_big_words)][1]).' '.$important_motifs[unpack('N',$pack_indice_big_words)][2]."\n";

		      push @{$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 0 ]}, $pack_indice_big_words;
		      push @{$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ]}, pack('n',$#{ $H_important_motifs{ $pack_indice_big_words }[ 1 ] });
		      
		      # ($not_verbose)or print "RECORD_TRINUC for $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp)." file ".$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_indice_big_words}[0]}[ $#{ $H_important_motifs{ $pack_indice_big_words }[ 1 ] } ])][0]."\n";

		      &RECORD_TRINUC( $#{ $H_important_motifs{ $pack_indice_big_words }[ 1 ] }, $num_w1, $num_w2, $r_pack_boxsp, \$pack_indice_big_words );
		      
		      # ($not_verbose)or print "nouvel enregistrement de mot deja existant $$r_chaine1 CAS2\ncompte vaut ".substr($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[ 1 ], $#{ $H_important_motifs{ $pack_indice_big_words }[ 1 ] }, 1)." FILE ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]."\n";
		      if(length($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[ 1 ])<$#{ $H_important_motifs{ $pack_indice_big_words }[ 1 ] }+1){ die "$prog_tag [Error] chaine trop courte, line ".__LINE__."\n"; }		    
		      
		      $bool_not_half_recorded = 0;
		      last;
		  }
	      }
	  }
	    
	    if($bool_not_half_recorded)
	    {
		# if we did not find the diff_pos, we have to record data related to trinuc1_trinuc2 and the counters
		# ($not_verbose)or print "WORD NOT PRESENT IN FIC_TRINUC BUT WE HAVE TO SEE IF IT IS PRESENT IN TRINUC1_TRINUC2_SPACER file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]."\n";

		my $bool_not_existing_H = -1;
		
		if($bool_not_fic_diff)
		{

		    for my $ind_H_m(0..$#{$H_important_motifs{ $pack_indice_big_words }[0]})
		    {
			if((unpack('n',${$H_important_motifs{ $pack_indice_big_words }[0]}[$ind_H_m]) == unpack('n',$$r_pack_file) )and( unpack('n',${$H_important_motifs{ $pack_indice_big_words }[1]}[$ind_H_m]) == unpack('n',$local_pack_pos_reelle_deb)))
			{
			    $bool_not_existing_H = $ind_H_m;
			    last;
			}
		    }
		}
		else
		{
		    $spacer_word1_word2_fic{ $pack_espace }{ $mot1 }{ $mot2 }{ $$r_pack_file } = $pack_indice_big_words;
		}

		push @{$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 0 ]}, $pack_indice_big_words;
		push @{$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[0]}, $pack_diff_pos;

		# print "on ajoute pour le fichier ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." pour les sous-mots $table_trinuc[$num_w1] ".unpack('U',$$r_pack_boxsp)." $table_trinuc[$num_w2] ln 1232 le mot ".$important_motifs[unpack('N',$pack_indice_big_words)][0].' '.unpack('U',$important_motifs[unpack('N',$pack_indice_big_words)][1]).' '.$important_motifs[unpack('N',$pack_indice_big_words)][2]."\n";

		if( $bool_not_existing_H == -1)
		{
		    push @{$H_important_motifs{ $pack_indice_big_words }[0]}, $$r_pack_file;
		    push @{$H_important_motifs{ $pack_indice_big_words }[1]}, $local_pack_pos_reelle_deb;
		    push @{$H_important_motifs{ $pack_indice_big_words }[2]}, $local_pack_pos_reelle_deb_autre_seq;
		    push @{$H_important_motifs{ $pack_indice_big_words }[3]}, $local_pack_shiftsp_autre_seq;
		    # ($not_verbose)or STDERR print "word already recorded, file ".$tab_ID_fic[unpack('n', $$r_pack_file)][0].' '.$important_motifs[unpack('N',$pack_indice_big_words)][0].' '.unpack('U',$important_motifs[unpack('N',$pack_indice_big_words)][0]).' '.$important_motifs[unpack('N',$pack_indice_big_words)][2].', we add positions pos1 '.unpack('n',$local_pack_pos_reelle_deb).' pos2 '.unpack('n',$local_pack_pos_reelle_deb).' shifts '.unpack('cc',$local_pack_shiftsp_autre_seq)." line ".__LINE__."\n";

		    push @{$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ]}, pack('n',$#{$H_important_motifs{ $pack_indice_big_words }[0]});

		    # ($not_verbose)or print "RECORD_TRINUC for $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp)." file ".$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_indice_big_words}[0]}[ $#{ $H_important_motifs{ $pack_indice_big_words }[ 1 ] } ])][0]."\n";

		    &RECORD_TRINUC( $#{ $H_important_motifs{ $pack_indice_big_words }[ 1 ] }, $num_w1, $num_w2, $r_pack_boxsp, \$pack_indice_big_words );
		    # ($not_verbose)or print "nouvel enregistrement de Himp et fictrinuc pour mot deja existant $$r_chaine1 file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." CAS31\n";

		}
		else
		{
		    push @{$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ]}, pack('n',$bool_not_existing_H);

		    # ($not_verbose)or print "RECORD_TRINUC for $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp)." file ".$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_indice_big_words}[0]}[ $bool_not_existing_H ])][0]."\n";

		    &RECORD_TRINUC( $bool_not_existing_H, $num_w1, $num_w2, $r_pack_boxsp, \$pack_indice_big_words );
		    # ($not_verbose)or print "nouvel enregistrement de fic_trinuc seulement pour mot deja existant $$r_chaine1 file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0]." CAS32\n";
		}

	    }
	    
	    if(! $not_verbose)
	    {
		my $verif_table_trinuc_num_w1 = $table_trinuc[$num_w1];
		my $verif_table_trinuc_num_w2 = $table_trinuc[$num_w2];
		&REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w1);
		&REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w2);
		if($important_motifs[ unpack('N',$pack_indice_big_words) ][0] !~ /$verif_table_trinuc_num_w1/){ die "$prog_tag [Error] dans remove w1 $mot1 sp ".unpack('U',$pack_espace)." w2 $mot2 ne contient pas r_w1 $table_trinuc[$num_w1] \nchaine1 $$r_chaine1 ID ".unpack('N',$pack_indice_big_words)."\nCONCATENATION $mot1 doit correspondre a ".($important_motifs[ unpack('N',$pack_indice_big_words) ][0])." file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0].", line ".__LINE__."\n"; } 
		if($important_motifs[ unpack('N',$pack_indice_big_words) ][2] !~ /$verif_table_trinuc_num_w2/){ die "$prog_tag [Error] dans remove w1 $mot1 sp ".unpack('U',$pack_espace)." w2 $mot2 ne contient pas r_w2 $table_trinuc[$num_w2] \nchaine1 $$r_chaine1 ID ".unpack('N',$pack_indice_big_words)."\nCONCATENATION $mot2 doit correspondre a ".($important_motifs[ unpack('N',$pack_indice_big_words) ][2])." file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0].", line ".__LINE__."\n"; } 
	    }
	}
    }

    if($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_indice_big_words }[1] == 0)
    {
	die "$prog_tag [Error] trinuc1_trinuc2_spacer  $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp).' '.unpack('N',$pack_indice_big_words)." 1 vide pour $important_motifs[ unpack('N',$pack_indice_big_words) ][0] ".unpack('U',$important_motifs[  unpack('N',$pack_indice_big_words) ][1])." $important_motifs[ unpack('N',$pack_indice_big_words) ][2], line ".__LINE__."\n";
    }

    # removing of old useless records in tables and hash tables when counter will reach 0
    # MAJ of hash table (dyads of trinucleotides), remove of those useless
    # pour tous les sous mots repertories (indice num, position...)

    # creation of the new access to the new dyad of words in the trinuc H table


    #    print "ind deb $local_ind_deb i $local_i\n";
    for my $i_supp($local_ind_deb..$local_i)
    {

	my $pack_id_big_table = $r_champs_num->[$i_supp];
	if( $i_supp != $bool_not_existing_diff_pos )	
	{
	   
	    # print "id_big_table ".unpack('N',$pack_id_big_table)."\n";
	    # print "pack_indice_big_words ".unpack('N',$pack_indice_big_words)."\n";
	    # print "champs pos reelle de $i_supp $r_champs_pos_reelle->[$i_supp]\n";
	    # print "local_pos_reelle $local_pos_reelle_deb\n";
	    # print "trinuc1trinuc2 w1 $table_trinuc[$num_w1] w2 $table_trinuc[$num_w2] sp ".unpack('U',$$r_pack_boxsp)." id_big_t ".unpack('N',$pack_id_big_table)." de posdanstrinuc $r_champs_ind_diff_pos->[$i_supp] vaut ".unpack('U',${ $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 0 ] }[ $r_champs_ind_diff_pos->[$i_supp] ])."\n";
	    # print "mot $important_motifs[ $pack_id_big_table ][0] ".unpack('U',$important_motifs[ unpack('N',$pack_id_big_table) ][1])." $important_motifs[ unpack('N',$id_big_table) ][2]\n";
	    
	    if( substr($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ], unpack('n',${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ]), 1) eq '' )
	    {
		die "$prog_tag [Error] position ".unpack('n',${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ])." de la chaine ".$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ]." sans rien ".substr($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ], unpack('n',${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ]) , 1)."\nconcerne trinuc de ".unpack('n',$$r_pack_file)." $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp)." num ".unpack('N',$pack_id_big_table)." indice  $r_champs_ind_diff_pos->[$i_supp], line ".__LINE__."\n";
	    }
	    
	    # ($not_verbose)or print "chaine vaut $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ]\n";

	    if( substr($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ], unpack('n',${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ]) , 1) == 0 )
	    {
		die "$prog_tag [Error] champs a decrementer egale a ".substr($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ], unpack('n',${$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ]) , 1)." (0) rr_trinuc1_trinuc2_spacer $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp).' '.unpack('N',$pack_id_big_table)." 1 idpos $r_champs_ind_diff_pos->[$i_supp] ] possouschaine ".unpack('n', ${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ]).", line ".__LINE__."\n";
	    }
	    
	    substr($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ], unpack('n',${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ]), 1)--;
	    
	    # ($not_verbose)or STDERR print "on decremente rr_trinuc1_trinuc2_spacer $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp).' '.unpack('N',$pack_id_big_table)." 1  ".unpack('n',${$fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ])." qui passe a $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ]: $important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][ 0 ] ".unpack('U',$important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][ 1 ])." $important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][ 2 ]\n";

	    if($trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ]  !~ /^\d+$/){ die "$prog_tag [Error] line ".__LINE__." negatif $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ]\n"; }
	    
	    if( $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ] == 0 )
	    {
		# ($not_verbose)or print "on SUPPRIME TOUT rr_trinuc1_trinuc2_spacer  $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp).' '.unpack('N',$pack_id_big_table)." car il n'y a plus d'enregistrement dans le tableau\n\n";

	      # --------------------------------------------------------------------------------------------
	      # changed 20140917
	      # line ...
	      # %{$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }} = ();

		# replaced by
		@{$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }} = ();
		# $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table } = undef;
	      # --------------------------------------------------------------------------------------------
	    }

	}
	else
	{
	    if((exists $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 0 ])and($#{ $trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 0 ] } == -1))
	    {
		%{$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }} = ();
		die "$prog_tag [Error] suppression a un endroit ou on n'est pas suppose en faire file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0].", line ".__LINE__."\n";
	    }

	    # OK
	    if(! $not_verbose)
	    {
		my $verif_table_trinuc_num_w2 = $table_trinuc[$num_w2];
		&REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w2);
		if($important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][2] !~ /$verif_table_trinuc_num_w2/)
		{ 
		    die "$prog_tag [Error] dans remove w1 $important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][0] sp ".unpack('U',$important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][1])." w2 $important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][2] ne contient pas r_w2 $table_trinuc[$num_w2] \nBOUCLE SUPP file ".$tab_ID_fic[unpack('n',$$r_pack_file)][0].", line ".__LINE__."\n"; 
		} 
	    }
	}
    }

    # if(($table_trinuc[$num_w1] eq 'ggc')and($table_trinuc[$num_w2] eq 'gcc')and(unpack('U',$$r_pack_boxsp) == 22))
    # {
#	print "\n\n";
#	for my $i_supp(reverse $local_ind_deb..$local_i)
#	{
#	    my $pack_id_big_table = $r_champs_num->[$i_supp];

#	    if(( unpack('N',$pack_id_big_table) != $bool_not_existing_diff_pos )and(exists ${$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ]}[ $r_champs_ind_diff_pos->[$i_supp] ]))	
#	    {
		
#		print "on a rr_trinuc1_trinuc2_spacer $table_trinuc[$num_w1] $table_trinuc[$num_w2] ".unpack('U',$$r_pack_boxsp).' '.unpack('N',$pack_id_big_table)." 1 de $r_champs_ind_diff_pos->[$i_supp] indposition ".unpack('n',${ $fic_trinuc1_trinuc2{ $$r_pack_file }{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }[ 1 ] }[ $r_champs_ind_fic_tri->[$i_supp] ])." qui est a ${$trinuc1_trinuc2_spacer{ $num_w1 }{ $num_w2 }{ $$r_pack_boxsp }{ $pack_id_big_table }[ 1 ]}[ $r_champs_ind_diff_pos->[$i_supp] ]: $important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][ 0 ] ".unpack('U',$important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][ 1 ])." $important_motifs[ unpack('N',$r_champs_num->[$i_supp]) ][ 2 ]\n";
#	    }
#	}
	# die;
    # }
}

# subprogramme dedicated to the extend of dyads when some of them overlap
sub EXTEND_DYAD(\$)
{
    my($r_ind_important_motifs) = @_;

    # ($not_verbose)or print "SUB EXTEND_DYAD @_\n";
    # for each file concerned
    for my $pack_file(keys %fic_trinuc1_trinuc2)
    {
	for my $num_w1(keys %{ $fic_trinuc1_trinuc2{$pack_file} })
	{
	    for my $num_w2(keys %{ $fic_trinuc1_trinuc2{$pack_file}{$num_w1} })
	    {
		for my $pack_boxsp(keys %{ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2} })
		{
		    my @champs_pack_num     = ();
		    my @champs_pos_reelle   = ();
		    my @champs_pos_reelle2  = ();
		    my @champs_ind_diff_pos = ();
		    my @champs_ind_fic_tri  = ();
		    my %champs_pos_doublet_trinuc = ();

		    if($#{ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }> 0)
		    {
			# print "***************************************************\n";
			# ($not_verbose)or STDERR print "COUPLE CONCERNE fic ".$tab_ID_fic[unpack('n',$pack_file)][0]." w1 $table_trinuc[$num_w1] w2 $table_trinuc[$num_w2] boxsp ".unpack('U',$pack_boxsp)."\n";

			my %H_redond = ();

			# pour chaque differente position dans H_important_motifs
			for my $cpt_w(0..$#{ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] } )
			{
#			    my @ref_to_increase = ();
#			    my $factor = 0;

			    # print "num ".unpack('N',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w]).' ';   
			    # print "postab ".unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[$cpt_w]).' ';   
			    # print "pos reelle ".unpack('n',${ $H_important_motifs{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w] }[1] }[ unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[$cpt_w]) ])."\n";

			    for my $several_diff_pos(0..$#{ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w] }[ 0 ] })
			    {
				
				my $unpack_n_fic_trinuc0 = unpack('N',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w]);
				my $unpack_n_fic_trinuc1 = unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[$cpt_w]);
				if((! exists $H_redond{"$unpack_n_fic_trinuc0 $unpack_n_fic_trinuc1 $several_diff_pos"})and
				   (defined  $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w] }[ 1 ] )and(length($trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w] }[ 1 ] )> $unpack_n_fic_trinuc1)and
				   (substr( $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w] }[ 1 ] , $unpack_n_fic_trinuc1, 1) != 0)
				   )
				{

				    $H_redond{"$unpack_n_fic_trinuc0 $unpack_n_fic_trinuc1 $several_diff_pos"} = 1;
			
				    # PROBLEME: il ne faut tenir compte pour la multiplication que des trinuc1_trinuc2 du meme fichier (donc dans file, necessite comptage preliminaire# ??!!!!!!!!!!
				    my $pack_num = ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w];
				    my $diff_pos_en_plus = unpack('U',${$trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$pack_num}[0]}[ $several_diff_pos ]);
				    my $pos_reelle = unpack('n',${ $H_important_motifs{ $pack_num }[1] }[ $unpack_n_fic_trinuc1 ]);

				    push @{ $champs_pos_doublet_trinuc{ $pos_reelle + $diff_pos_en_plus }[0] }, $pack_num; # @champs_num

				    push @{ $champs_pos_doublet_trinuc{ $pos_reelle + $diff_pos_en_plus }[1] }, unpack('n',${ $H_important_motifs{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w] }[1] }[ $unpack_n_fic_trinuc1 ]); # @champs_pos_reelle seq 1

				    push @{ $champs_pos_doublet_trinuc{ $pos_reelle + $diff_pos_en_plus }[2] }, unpack('n',${ $H_important_motifs{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w] }[2] }[ $unpack_n_fic_trinuc1 ]); # @champs_pos_reelle seq 2

				    push @{ $champs_pos_doublet_trinuc{ $pos_reelle + $diff_pos_en_plus }[3] }, $several_diff_pos; # @champs_ind_diff_pos

				    push @{ $champs_pos_doublet_trinuc{ $pos_reelle + $diff_pos_en_plus }[4] }, $cpt_w; # @champs_ind_fic_tri

				    # record num ID of the dyad word in the global table important_motifs 

				    # print "num $unpack_n_fic_trinuc0 ";   
				    # print "postab $unpack_n_fic_trinuc1 ";   
				    # print "pos reelle ".unpack('n',${ $H_important_motifs{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[$cpt_w] }[1] }[ $unpack_n_fic_trinuc1 ]);
				    # print " inddiffpos $several_diff_pos ";
				    # print "ind_fic_tri $cpt_w\n";
				}
			    }
			}
		    }
		    
		    if(scalar(keys %champs_pos_doublet_trinuc )!= 0)
		    {

			# we first sort by trinucleotid doublet position (each corresponding to 4 tables (0, 1, 2, 3) to sort)
			for my $cle_cpdt(sort numerique keys %champs_pos_doublet_trinuc)
			{

			    # tri en fonction de la position
			    # &tri_rapide(\@champs_pos_reelle, \@champs_num, \@champs_ind_diff_pos, \@champs_ind_fic_tri, 0, scalar(@champs_num)-1);
			
			    &tri_rapide(  0, scalar(@{ $champs_pos_doublet_trinuc{$cle_cpdt}[1] })-1, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[1] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[0] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[2] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[3] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[4] });

			    my $deb_tri = 0;
			    my $fin_tri = 0;
			    for my $interv_tri(0..$#{ $champs_pos_doublet_trinuc{$cle_cpdt}[1] }-1)
			    {
				if( ${ $champs_pos_doublet_trinuc{$cle_cpdt}[1] }[$interv_tri] != ${ $champs_pos_doublet_trinuc{$cle_cpdt}[1] }[$interv_tri+1] )
				{
				    $deb_tri = $interv_tri + 1;
				    ($deb_tri == $fin_tri)or &tri_rapide(  $deb_tri, $fin_tri, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[2] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[0] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[1] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[3] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[4] });
				}
				else
				{
				    # $fin_tri = $interv_tri + 1; CHANGED 010506
				    $fin_tri = $interv_tri;
				    ($fin_tri != $#{ $champs_pos_doublet_trinuc{$cle_cpdt}[1] })or($deb_tri == $fin_tri)or &tri_rapide(  $deb_tri, $fin_tri, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[2] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[0] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[1] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[3] }, \@{ $champs_pos_doublet_trinuc{$cle_cpdt}[4] });
				}
			    }

			    # ??AMELIORABLE, RECOPIE INUTILE?
			    push @champs_pack_num     , @{ $champs_pos_doublet_trinuc{$cle_cpdt}[0] }; 
			    push @champs_pos_reelle   , @{ $champs_pos_doublet_trinuc{$cle_cpdt}[1] }; 
			    push @champs_pos_reelle2  , @{ $champs_pos_doublet_trinuc{$cle_cpdt}[2] }; 
			    push @champs_ind_diff_pos , @{ $champs_pos_doublet_trinuc{$cle_cpdt}[3] }; 
			    push @champs_ind_fic_tri  , @{ $champs_pos_doublet_trinuc{$cle_cpdt}[4] }; 
			}

			# display to check tables

			# if(! $not_verbose)
			# {
			  #  print "DEB TAB\n";
			  #  for my $i(0..$#champs_pack_num)
			  #  {
			#	print "champsnum ".unpack('N',$champs_pack_num[$i])." mot1 $important_motifs[unpack('N',$champs_pack_num[$i])][0] mot2 $important_motifs[ unpack('N',$champs_pack_num[$i]) ][2], position1 $champs_pos_reelle[$i], position2 $champs_pos_reelle2[$i],";
			#	print " ind fictri $champs_ind_fic_tri[$i] ";
			#	print " ind diff_pos $champs_ind_diff_pos[$i] ";
			#	print ", la position dans le H_impo vaut ".unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[$champs_ind_fic_tri[$i]]).", ";
			#	print "doit donner la position ";

			#	print unpack('n',${ $H_important_motifs{$champs_pack_num[$i]}[1] }[ unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[$champs_ind_fic_tri[$i]]) ])." autre position ".unpack('n',${ $H_important_motifs{$champs_pack_num[$i]}[2] }[ unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[$champs_ind_fic_tri[$i]]) ])."\n";
			#        print "nboccu dans trinuc1_trinuc2_spacer ".substr($trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i]}[1], unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[$champs_ind_fic_tri[$i]]), 1)."\n"; 
			#    }
			# }
		    }

		    %champs_pos_doublet_trinuc = ();
		
		    # NE PAS REMETTRE, FAUX
		    # if($#champs_pack_num > 1)
		    
		    # comparison of the different position so as to fusion dyads if necessary
		
		    my $chaine1 = '';
		    my $chaine2 = '';
		    
		    my $pack_pos_reelle_deb           = ();
		    my $pack_pos_reelle_deb_autre_seq = ();
		    my $pack_shiftsp_autre_seq        = ();
		    my $ind_deb                  = -1;
		    my $cumul_diff_dep;
		    
		    # print "file ".$tab_ID_fic[unpack('n',$pack_file)][0]."\n";

		    for my $i(0..$#champs_pack_num-1){
			
			my $champs_num = unpack('N',$champs_pack_num[$i]);
			my $champs_num_plus1 = unpack('N',$champs_pack_num[$i+1]);
			# ($not_verbose)or print "champs de $i: num $champs_num pos1 $champs_pos_reelle[$i] pos2 $champs_pos_reelle2[$i]\n";
			
			if(! $not_verbose)
			{
			    my $verif_table_trinuc_num_w1 = $table_trinuc[$num_w1];
			    my $verif_table_trinuc_num_w2 = $table_trinuc[$num_w2];
			    &REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w1);
			    &REMPLACE_N_PAR_ANTISLASHW(\$verif_table_trinuc_num_w2);
			    if($important_motifs[ $champs_num ][0] !~ /$verif_table_trinuc_num_w1/){ die "$prog_tag [Error] num $champs_num\n$important_motifs[ $champs_num ][0] ".unpack('U',$important_motifs[ $champs_num ][1])." $important_motifs[ $champs_num ][2] ne contient pas w1 $table_trinuc[$num_w1] dans extenddyad apres tri, line ".__LINE__."\n";}			    
			    if($important_motifs[ $champs_num_plus1 ][0] !~ /$verif_table_trinuc_num_w1/){ die "$prog_tag [Error] num $champs_num_plus1\n$important_motifs[ $champs_num_plus1 ][0]  ".unpack('U',$important_motifs[ $champs_num_plus1 ][1])." $important_motifs[ $champs_num_plus1 ][2] ne contient pas w1 $table_trinuc[$num_w1] dans extenddyad apres tri, line ".__LINE__."\n";}			    
			    if($important_motifs[ $champs_num ][2] !~ /$verif_table_trinuc_num_w2/){ die "$prog_tag [Error] num $champs_num\n$important_motifs[ $champs_num ][0]  ".unpack('U',$important_motifs[ $champs_num ][1])." $important_motifs[ $champs_num ][2] ne contient pas w2 $table_trinuc[$num_w2] dans extenddyad apres tri, line ".__LINE__."\n";}			    
			    if($important_motifs[ $champs_num_plus1 ][2] !~ /$verif_table_trinuc_num_w2/){ die "$prog_tag [Error] num $champs_num_plus1\n$important_motifs[ $champs_num_plus1 ][0]  ".unpack('U',$important_motifs[ $champs_num_plus1 ][1])." $important_motifs[ $champs_num_plus1 ][2] ne contient pas w2 $table_trinuc[$num_w2] dans extenddyad apres tri, line ".__LINE__."\n";}			    
			}

			my $diff_dep  = $champs_pos_reelle[$i+1]  - $champs_pos_reelle[$i];

			defined ${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i]}[0] }[ $champs_ind_diff_pos[$i] ] or die "$prog_tag [Error] trinuc1_trinuc2_spacer of num_w1:$num_w1 (corresponds to $table_trinuc[$num_w1]) of num_w2:$num_w2 (corresponds to $table_trinuc[$num_w2]) of pack_boxsp:$pack_boxsp of (champs_pack_num de i:$i :$champs_pack_num[$i]) of 0 (array_index) of (champs_ind_diff_pos de i:$i :$champs_ind_diff_pos[$i] (array index) ) not defined line ".__LINE__."\n";

			# my $diff_dep2 = $champs_pos_reelle2[$i+1] - $champs_pos_reelle2[$i];
			my $diff_diff_pos = ( $champs_pos_reelle[$i+1] 
					      + unpack('U',${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i+1]}[0] }[ $champs_ind_diff_pos[$i+1] ]) 
					      - $champs_pos_reelle[$i] 
					      - unpack('U',${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i]}[0] }[ $champs_ind_diff_pos[$i] ]) );
			my $diff_diff_pos2 = ( $champs_pos_reelle2[$i+1] 
					      + unpack('U',${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i+1]}[0] }[ $champs_ind_diff_pos[$i+1] ]) 
					      - $champs_pos_reelle2[$i] 
					      - unpack('U',${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i]}[0] }[ $champs_ind_diff_pos[$i] ]) );

			my ($shitsp11,$shiftsp10) = unpack('cc',${ $H_important_motifs{ $champs_pack_num[$i+1] }[3] }[ unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[ $champs_ind_fic_tri[$i+1] ] ) ]);
			my $diff_shift = $shiftsp10 - $shitsp11;
			($shitsp11,$shiftsp10) = unpack('cc',${ $H_important_motifs{ $champs_pack_num[$i] }[3] }[ unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[ $champs_ind_fic_tri[$i] ] ) ]);
			$diff_shift -= ($shiftsp10 - $shitsp11);
			# print "diff_pos $diff_diff_pos resulte de champsposreelleiplus1 $champs_pos_reelle[$i+1] plus unpackUtrinuc1trinuc2spaceriplus1 ".unpack('U',${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i+1]}[0] }[ $champs_ind_diff_pos[$i+1] ])." moins champposreellei $champs_pos_reelle[$i] moins unpackUde trinuc1trinuc2spaceri ".unpack('U',${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i]}[0] }[ $champs_ind_diff_pos[$i] ])."\n";
			# print "diff_pos2 $diff_diff_pos2 resulte de champsposreelle2iplus1 $champs_pos_reelle2[$i+1] plus unpackUtrinuc1trinuc2spaceriplus1 ".unpack('U',${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i+1]}[0] }[ $champs_ind_diff_pos[$i+1] ])." moins champposreelle2i $champs_pos_reelle2[$i] moins unpackUde trinuc1trinuc2spaceri ".unpack('U',${ $trinuc1_trinuc2_spacer{$num_w1}{$num_w2}{$pack_boxsp}{$champs_pack_num[$i]}[0] }[ $champs_ind_diff_pos[$i] ])."\n";
			# we get beginning of the first word of new dyad
#			if(($chaine1 eq '')and((my $diff_w1 = length($important_motifs[ $champs_pack_num[$i] ][0]) - $diff_dep ) >= 0)and($diff_diff_pos == 0))

			if(($chaine1 eq '')and($diff_diff_pos == 0)and($diff_diff_pos2 == 0)and($diff_shift == 0))
			{
			    my $diff_w1 = length($important_motifs[ unpack('N',$champs_pack_num[$i]) ][0]) - $diff_dep;
			    # ($not_verbose)or print "CAS 0\n";
			    $ind_deb = $i;
			    my $spa = '-' x unpack('U',$important_motifs[ unpack('N',$champs_pack_num[$i]) ][1]);
			    $chaine1 = $important_motifs[ unpack('N',$champs_pack_num[$i]) ][0].$spa.$important_motifs[ unpack('N',$champs_pack_num[$i]) ][2];

			    $spa = '-' x unpack('U',$important_motifs[ unpack('N',$champs_pack_num[$i+1]) ][1]);
			    my $tiret_dep = '-' x $diff_dep;
			    $chaine2 = $tiret_dep.$important_motifs[ unpack('N',$champs_pack_num[$i+1]) ][0].$spa.$important_motifs[ unpack('N',$champs_pack_num[$i+1]) ][2];
			    if( (my $diff_w = length($chaine2) - length($chaine1) ) > 0 )
			    { 
				$chaine1 .= '-'x $diff_w; 
			    }
			    elsif($diff_w < 0)
			    {
				$chaine2 .= '-'x -$diff_w;
			    }
			    
			    $cumul_diff_dep = $diff_dep;
			    
			    # ($not_verbose)or print "chaine1 $chaine1 pour ".unpack('N',$champs_pack_num[$i])."\n";
			    # ($not_verbose)or print "chaine2 $chaine2 pour ".unpack('N',$champs_pack_num[$i+1])."\n";
			    
			    for my $lettre(0..length($chaine1)-1)
			    {
				(substr($chaine1,$lettre,1) ne '-')or substr($chaine1,$lettre,1) = substr($chaine2,$lettre,1);
				if( (substr($chaine2,$lettre,1) ne '-')and(substr($chaine2,$lettre,1) ne substr($chaine1,$lettre,1)) )
				{
				    die "$prog_tag [Error] chaine1 $chaine1 et chaine2 $chaine2 ne correspondent pas!!!\n".substr($chaine2,$lettre,1)." ne ".substr($chaine1,$lettre,1)." lettre $lettre, line ".__LINE__."\n";
				}
			    }
			    
			    # position of the motif we will record
			    $pack_pos_reelle_deb = pack('n',$champs_pos_reelle[$i]);
			    $pack_pos_reelle_deb_autre_seq = pack('n',$champs_pos_reelle2[$i]);
#			    $pack_pos_reelle_deb_autre_seq = ${ $H_important_motifs{ $champs_pack_num[$i] }[2] }[ unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[ $champs_ind_fic_tri[$i] ]) ];

			    # we don't need to unpack as values are not evaluated, we will only need to copy the value
			    $pack_shiftsp_autre_seq  = ${ $H_important_motifs{ ${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[0] }[ $champs_ind_fic_tri[$i] ] }[3] }[ unpack('n',${ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2}{$pack_boxsp}[1] }[ $champs_ind_fic_tri[$i] ]) ];

			    if($i+1 == $#champs_pack_num)
			    {
#				if(($chaine1 eq "cgggaa-----------------tggttg")and($pack_file eq "SCO2634.Rv2466c")){ die "bien trouve"; } 
				# ($not_verbose)or print "enregistrement CAS 0\n";
				# ($not_verbose)or print "chaine1 fixee $chaine1\n";

				&REMOVE($pack_pos_reelle_deb, $pack_pos_reelle_deb_autre_seq, $pack_shiftsp_autre_seq, $ind_deb,$i+1,$diff_w1,$num_w1,$num_w2,\$pack_file,\$pack_boxsp,\@champs_pack_num,\@champs_ind_diff_pos,\@champs_ind_fic_tri, \@champs_pos_reelle, \$chaine1,$r_ind_important_motifs);
				# initialization of variable for other words
				$chaine1 = '';
				$chaine2 = '';
				$pack_pos_reelle_deb = ();
				$pack_pos_reelle_deb_autre_seq = ();
			    }
			}
			elsif(($diff_diff_pos == 0)and($diff_diff_pos2 == 0)and($diff_shift == 0)) 
			{
			    my $diff_w1 = length($important_motifs[ unpack('N',$champs_pack_num[$i]) ][0]) - $diff_dep;

			    my $spa = '-' x unpack('U',$important_motifs[unpack('N', $champs_pack_num[$i+1]) ][1]);
			    
			    # ($not_verbose)or print "CAS 1\n";

			    $cumul_diff_dep += $diff_dep;
			    my $tiret_dep = '-' x $cumul_diff_dep;
			    $chaine2 = $tiret_dep.$important_motifs[ unpack('N',$champs_pack_num[$i+1]) ][0].$spa.$important_motifs[ unpack('N',$champs_pack_num[$i+1]) ][2];
			    if((my $diff_w = length($chaine2) - length($chaine1)) > 0)
			    { 
				$chaine1 .= '-'x $diff_w;
			    }
			    elsif($diff_w < 0)
			    {
				$chaine2 .= '-'x -$diff_w;
			    }
			    # ($not_verbose)or print "chaine1 $chaine1\n";
			    # ($not_verbose)or print "chaine2 $chaine2 pour ".unpack('N',$champs_pack_num[$i+1])."\n";
			
			    for my $lettre(0..length($chaine1)-1)
			    {
				(substr($chaine1,$lettre,1) ne '-')or substr($chaine1,$lettre,1) = substr($chaine2,$lettre,1);
			    }
			    
			    if($i+1 == $#champs_pack_num)
			    { 
				# ($not_verbose)or print "enregistrement CAS 1\n";
				# ($not_verbose)or print "chaine1 fixee $chaine1\n";

				&REMOVE($pack_pos_reelle_deb, $pack_pos_reelle_deb_autre_seq, $pack_shiftsp_autre_seq, $ind_deb,$i+1,$diff_w1,$num_w1,$num_w2,\$pack_file,\$pack_boxsp,\@champs_pack_num,\@champs_ind_diff_pos,\@champs_ind_fic_tri, \@champs_pos_reelle, \$chaine1,$r_ind_important_motifs);
				# initialization of variable for other words
				$chaine1 = '';
				$chaine2 = '';
				$pack_pos_reelle_deb = ();
				$pack_pos_reelle_deb_autre_seq = ();
			    }
			}
			elsif($chaine1 ne '')
			{

			    my $diff_w1 = length($important_motifs[ unpack('N',$champs_pack_num[$i]) ][0]) - $diff_dep;
			    # ($not_verbose)or print "CAS 2 (enregistrement)\n";
			    # ($not_verbose)or print "chaine1 fixee $chaine1\n";

			    # print "diffw1 pas bon diffw1 vaut $diff_w1 car\nlength($important_motifs[ unpack('N',$champs_pack_num[$i]) ][0]) - ($champs_pos_reelle[$i+1] - $champs_pos_reelle[$i]) inf a 0\n";

			    # ($not_verbose)or print "chaine2 $chaine2\n";

			    &REMOVE($pack_pos_reelle_deb, $pack_pos_reelle_deb_autre_seq, $pack_shiftsp_autre_seq, $ind_deb,$i,$diff_w1,$num_w1,$num_w2,\$pack_file,\$pack_boxsp,\@champs_pack_num,\@champs_ind_diff_pos,\@champs_ind_fic_tri, \@champs_pos_reelle, \$chaine1,$r_ind_important_motifs);
			    
			    # initialization of variable for other words
			    $chaine1 = '';
			    $chaine2 = '';
			    $pack_pos_reelle_deb = ();
			    $pack_pos_reelle_deb_autre_seq = ();
			}
			
			else
			{

			    # ($not_verbose)or print "CAS 3 (reinitialisation des variables)\n";
			    # print "diffw1 pas bon diffw1 vaut $diff_w1 car\nlength($important_motifs[ unpack('N',$champs_pack_num[$i]) ][0]) - ($champs_pos_reelle[$i+1] - $champs_pos_reelle[$i]) inf a 0\n";
			    # die "diffw1 pas bon diffw1 vaut $diff_w1 car\nlength($important_motifs[ unpack('N',$champs_pack_num[$i]) ][0]) - ($champs_pos_reelle[$i+1] - $champs_pos_reelle[$i]) inf a 0\n";
			    
			    $chaine1 = '';
			    $chaine2 = '';
			    $pack_pos_reelle_deb = ();
			    $pack_pos_reelle_deb_autre_seq = ();
			    # print "\n$important_motifs[ unpack('N',$champs_pack_num[$i]) ][0] ".unpack('U',$important_motifs[ unpack('N',$champs_pack_num[$i]) ][1])." $important_motifs[ unpack('N',$champs_pack_num[$i]) ][2]\net\n$important_motifs[ unpack('N',$champs_pack_num[$i+1]) ][0] ".unpack('U',$important_motifs[ unpack('N',$champs_pack_num[$i+1]) ][1])." $important_motifs[ unpack('N',$champs_pack_num[$i+1]) ][2]\npas assez proches\nlongueur ".length($important_motifs[ unpack('N',$champs_pack_num[0]) ][0])." pos reelle 1 ".$champs_pos_reelle[$i]." pos relle 2 ".$champs_pos_reelle[$i+1]."\n";
			}
		    
		    }
		
		}
		%{ $fic_trinuc1_trinuc2{$pack_file}{$num_w1}{$num_w2} } = ();  
	    }    
	}
    }
}

# ***************************************************************************************************************** #

# *****************************************************************************************************************
# subprogramme to cross results obtained for each ortholog pair
# param 1: ref to H table which allow to recover a motif thanks according the trinucleotid dyads it is composed by
# param 2: ref to table of every interesting motifs (-35 box spacer and -10 box)
# param 3: ref to H table associated to prec table: foreach motif are recorded File -35 position, spacershift for
#   first seq and spacershift for second seq (other sp)
# param 4: ref to table where are recorded allowed spacershifts

# GLOBAL PARAM: not_verbose, infos
sub RES_COMMON_ORTHO
{
    $prev_time_long_or_not = time();

    # ($not_verbose)or print "SUB RES_COMMON_ORTHO\n";
    if($bool_print_cross)
    {
		open(FCROSS,"> CROSS$name_fic_sp[0]"."_$name_fic_sp[1]".'.txt')or die "$prog_tag [Error] Impossible to create CROSS$name_fic_sp[0]"."_$name_fic_sp[1]".'.txt'.":$!, line ".__LINE__."\n";
		print FCROSS "$infos\n";
	# ($not_verbose)or print "\n\nCROISEMENT\n"; # "
    }
    if($bool_print_sp_align){

	open(FALIGNSP1,"> $align_id1")or die "$prog_tag [Error] Impossible to create $align_id1 :$!, line ".__LINE__."\n";
	open(FALIGNSP2,"> $align_id2")or die "$prog_tag [Error] Impossible to create $align_id2 :$!, line ".__LINE__."\n";
	print FALIGNSP1 $infos;
	print FALIGNSP2 $infos;

	# ($not_verbose)or print "ALIGNEMENT\n";
    }



    my $indice_subsp = int($#spacer_shift_fic / 2); # indice du milieu de tableau

    # TRAITEMENT DE LOTS DE TRINUCS, PAS DE TOUT ****************************


    # tri par ordre decroissant du nombre de mots concernes...
    %H_nbw = ();
    
    while( my ($num_subw1) = each %trinuc1_trinuc2_spacer){

		# ($not_verbose)or print "subw1 $rev_H_trinuc[$num_subw1] parcours r_trinuc1_trinuc2_spacer line ".__LINE__."\n";
		
		while( my ($num_subw2) = each %{$trinuc1_trinuc2_spacer{$num_subw1}}){
	
		    &TEST_NBSEQ_PER_INTERV($trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}, $indice_subsp, \%H_nbw, $num_subw1, $num_subw2);
		}
    }
    
    print "H_nbw recorded\n";
    # in the order of the number of occurrences for the dyad of trinucleotides

    # THIS REVERSE is IMPORTANT
    # sort from the motif with the highest number of occurrences to this with the lowest

    for my $nbw(reverse sort numerique keys %H_nbw){

	if($bool_print_cross){
 
	    print FCROSS "$nbw occurrences\n"; 
	    #  ($not_verbose)or print "$nbw occurrences\n"; 
	}

	@seq_sp               = (); # table of parts of sp sequences
	@ID_f                 = (); # table of file IDs
	@infos_div            = (); # record data needed by get2.. subprog to avoid redondance in seq of each sp
	%prev_infos_div       = (); # remember data useded by get2.. subprog to avoid redondance in seq of each sp
	$entete               = ''; # string before results when we write in species files
	$intitules_ind_prem_l = ''; # record for the concerned line (number= key) the title
	
	@pos_tri              = (); # record pos of trinuc -35 sp1 -10 sp1 -35 sp2 -10 sp2


	# for my $verif_H(keys  %{$H_nbw{$nbw}})
	    # {
	      #  print 'verifH ';
	      #  print $$verif_H;
	      #  print $table_trinuc[$verif_H], "\n";
	    # }

	    for my $num_subw1( map { $_->[1] }
			       reverse sort { $a->[0] <=> $b->[0] }
			       map { [ length($table_trinuc[$_]), $_ ] } keys %{$H_nbw{$nbw}})
	    {
		my $subw1 = $table_trinuc[$num_subw1];
		$LOCAL_WIDTH_OF_L_HBOXES = length($subw1);

		# in the reverse order of the length of the second word 
		for my $num_subw2(map { $_->[1] }
				  reverse sort { $a->[0] <=> $b->[0] }
				  map { [ length($table_trinuc[$_]), $_ ] } keys %{$H_nbw{$nbw}{$num_subw1}})
		{
		    my $subw2 = $table_trinuc[$num_subw2];
		    $LOCAL_WIDTH_OF_R_HBOXES = length($subw2);

		    # in the order from the smaller to the higher spacer between the two boxes

		    # boolean to know if we have to set a new line between results
		    my $bool_saut_ln = 0;

		    @seq_sp               = (); # contains sequence part of sp1
		    @ID_f                 = (); # stocks IDs given by GET2.. subprog
		    @infos_div            = (); # contains ( ID sp1, pos in sp1, ID sp2, pos in sp2)
		    %prev_infos_div       = (); # contains (keys "last ID sp1 last pos in sp1" or "last ID sp2 last pos in sp2")
		    $entete               = '';
		    $intitules_ind_prem_l = "cas subw1 $subw1 subw2 $subw2 subsp";
		    @pos_tri              = ();

		    for my $pack_subsp(sort unpack_numerique_U keys %{$H_nbw{$nbw}{$num_subw1}{$num_subw2}})
		    {
			# display with shift in spacer between files

			# we use trinuc1_trinuc2_spacer and not r_H_nbw because we group related results (some are not present in r_H_nbw as they are not interesting like this, they are interesting as they are close to interesting results)
			if(scalar(keys( %{$trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}})) > 0)
			{
			    if($bool_print_cross)
			    {
				print FCROSS "cas subw1 $subw1 subw2 $subw2 subsp ".unpack('U',$pack_subsp)."\n";
			    }
			    $intitules_ind_prem_l .= ' '.unpack('U', $pack_subsp);
			    ($not_verbose)or print "cas subw1 $subw1 subw2 $subw2 subsp ".unpack('U',$pack_subsp)."\n";
			}
			else
			{ 
			    ($not_verbose)or print "we do not have enough values in r_trinuc1_trinuc2_spacer line ".__LINE__." for subw1 $subw1 subw2 $subw2 subsp ".unpack('U',$pack_subsp)."\n";
			    next; 
			}

			for my $pack_num(reverse sort unpack_numerique_N keys %{$trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}})
			{
			    if($#{$trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}{$pack_num}[0]} != -1)
			    {
				for my $id_pos1(0..$#{$trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}{$pack_num}[0]})
				{

				    (defined $trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}{$pack_num}[1])or next;

				    my $pos1 = unpack('U',${$trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}{$pack_num}[0]}[ $id_pos1 ]);
				    if(defined $pos1)
				    {
					my $line = '';
					my $num = unpack('N',$pack_num);
					my $spa = '-'x unpack('U',$important_motifs[ $num ][1]);
					my $pos2 = $pos1 + length($subw1) + unpack('U',$pack_subsp);

					$line = $important_motifs[$num][0].$spa.$important_motifs[ $num ][2];

					if(! $not_verbose)
					{
					    my $subw1_nByw = $subw1;
					    &REMPLACE_N_PAR_ANTISLASHW(\$subw1_nByw);
					    if(substr($line, $pos1, length($subw1)) !~ /$subw1_nByw/)
					    {
						die "$prog_tag [Error] PB line ".__LINE__.', '.(substr($line, $pos1, length($subw1)))." ne $subw1\n";
					    }
					}

					if($bool_print_cross)
					{
					    # print "on veut pos1 $pos1 pos2 $pos2 dans\nline ".__LINE__."\n";
					    substr($line, $pos1, length($subw1)) = uc($subw1);
					    substr($line, $pos2, length($subw2)) = uc($subw2);
					}

					for my $IDfic(0..$#{$H_important_motifs{$pack_num}[0]})
					{
					    if( (length($trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}{$pack_num}[1])> $IDfic)and(substr( $trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}{$pack_num}[1], $IDfic, 1) ne '')and( substr( $trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}{$pack_num}[1], $IDfic, 1) != 0 ) )
					    {
						$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][0] =~ /^(.+?)__(.+?)$/;

						my $unpack_H_imp_1 = unpack('n',${$H_important_motifs{$pack_num}[1]}[$IDfic]);
						my $unpack_H_imp_2 = unpack('n',${$H_important_motifs{$pack_num}[2]}[$IDfic]);
					
						# @infos_div = ($2, $unpack_H_imp_2+$pos1, $1, $unpack_H_imp_1+$pos1); 

						# print "dol1 $1\n";
						# print "unpack_H_imp_2 $unpack_H_imp_2\n";
						# print "pos1 $pos1\n";
						# print 'r_tab_ID_fic '.$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][1]." -1\n";
						# print "shift_between_end_seq_and_deb_trad $shift_between_end_seq_and_deb_trad\n";
						# print "dol2 $2\n";
						# print "unpack_H_imp_1 $unpack_H_imp_1\n";
						# print 'r_tab_ID_fic '.$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][2]." -1\n";

						# SI MEME ORDRE POUR ID DE FICHIERS ET SEQUENCES QU'IL CONTIENT
						
						@infos_div = ($1, $unpack_H_imp_2 + $pos1 - $tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][1] - 1 + $shift_between_end_seq_and_deb_trad, $2,  $unpack_H_imp_1 + $pos1 - $tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][2] - 1 + $shift_between_end_seq_and_deb_trad);
						
						# SI ORDRE DIFF POUR ID DE FICHIERS ET SEQUENCES QU'IL CONTIENT
						# @infos_div = ($2, $unpack_H_imp_2 + $pos1 - $tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][1] - 1 + $shift_between_end_seq_and_deb_trad, $1,  $unpack_H_imp_1 + $pos1 - $tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][2] - 1 + $shift_between_end_seq_and_deb_trad);
						
						$line .= " $tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][0] pos1 ".($unpack_H_imp_2 - $tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][1]  -1 + $shift_between_end_seq_and_deb_trad).' pos2 '.($unpack_H_imp_1 - $tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][2] - 1 + $shift_between_end_seq_and_deb_trad);
						
						# print "file $tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][0] subw1 $subw1 subw2 $subw2 subsp ".unpack('U',$pack_subsp)."\n";
						# we GET SEQ FOR EACH SPECIES to align them according to their similarities
						# ( ref to table of sp1 seq, motif position in sp1, spacershift + length to extract for sp1, ref to table of sp2 seq, motif position in sp2, spacershift + length to extract for sp2)
						my ($shift1,$shift2) = unpack('cc',${$H_important_motifs{$pack_num}[3]}[$IDfic]); # =~ /^(.+?) (.+?)$/;
						
						($not_verbose)or print "motif ".$important_motifs[ $num ][0].(unpack('U',$important_motifs[ $num ][1])).$important_motifs[ $num ][2].' file '.$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][0]."\n";
						# die "shift1 $shift1 shift2 $shift2\n";
						# ($not_verbose)or print "b_min2 vaut unpack_H_imp_1 $unpack_H_imp_1 - $LEFT_RIGHT_BOARDS_WIDTHS + pos $pos1 - 1 egale ".($unpack_H_imp_1 - $LEFT_RIGHT_BOARDS_WIDTHS + $pos1 - 1)."\n";
						
						
						($not_verbose)or print "traitement du motif $num correspondant au motif ".$important_motifs[ $num ][0].(unpack('U',$important_motifs[ $num ][1])).$important_motifs[ $num ][2].' file '.$tab_ID_fic[ unpack('n',${$H_important_motifs{$pack_num}[0]}[$IDfic]) ][0]." IDfic $IDfic\n";

						my $APPROX_SEQ_LENGTH_FOR_STAT = 2 * $LEFT_RIGHT_BOARDS_WIDTHS + $LOCAL_WIDTH_OF_L_HBOXES + $LOCAL_WIDTH_OF_R_HBOXES;

						# print "\nDEB GET2SEQFASTAMINMAX ******************************\n";
						&GET2SEQFASTAMINMAX(\@{$H_important_motifs{$pack_num}}, $IDfic, $important_motifs[ $num ], $unpack_H_imp_2 - $LEFT_RIGHT_BOARDS_WIDTHS + $pos1 - 1, $APPROX_SEQ_LENGTH_FOR_STAT + unpack('U',$pack_subsp) + $shift1 - $shift2, $unpack_H_imp_1 - $LEFT_RIGHT_BOARDS_WIDTHS + $pos1 - 1, $APPROX_SEQ_LENGTH_FOR_STAT + unpack('U',$pack_subsp), $LEFT_RIGHT_BOARDS_WIDTHS, \@infos_div , \%prev_infos_div );
						# print "FIN GET2SEQFASTAMINMAX ******************************\n\n";

					    }
					}

					if($bool_print_cross)
					{
					    $spa = ' 'x ($ordre_maxi+1-$pos1);
					    $line = $spa.$line;
					    substr($line, 16, 2) = unpack('U',$important_motifs[$num][1]);

					    print FCROSS "$line\n";
					    # ($not_verbose)or print $line."\n";
					}
				    }
				}
			    }
			}
			if( scalar(keys( %{$trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}{$pack_subsp}} )) > 0 )
			{
			    $bool_saut_ln = 1;
			}

		    } # end for my $pack_subsp
		    if(($bool_print_cross)and( $bool_saut_ln ))
		    {
				print FCROSS "\n";
		    }


		    if(! $not_verbose){
		      
			    print "verif seq\n";
				for my $i(0..1){
			
					print "espece $i\n";
					for(@{$seq_sp[$i]}){
					  
					  print "$_\n";
					}
			    }
		      
				# verif pos_tri
				for my $i(0..$#pos_tri){
			
					for my $j(0..$#{ $pos_tri[$i] }){
					  
					  print "postri $i $j vaut @{$pos_tri[$i][$j]}\n";
					}
			    }
		      print "end verif seq\n";
		      
		    }
		    # if(($subw1 eq 'gga')and($subw2 eq 'gtt'))
		    # {
			# die "FIN affichage seq correspondant a sigR\n";
		    # }

		    # **************************************************************************************
		    # TRI PAR LETTRE
		    @pos_tri_tmp = ();

		    # pour recuperation des lettres
		    @tab_lettre = ();

		    #  used to know lower and higher limit for spacer between trinuc boxes

		    # initialization of position of letters compared
		    for $sp(0..1){

			for my $i(0..$#{$pos_tri[$sp]})
			{ 
			    # shift is only to compense the first get of each letter
			    $pos_tri_tmp[$sp][$i][0] = $pos_tri[$sp][$i][0];   # pos at left of -35 box #**_____***
			    $pos_tri_tmp[$sp][$i][1] = $pos_tri[$sp][$i][0]+ $LOCAL_WIDTH_OF_L_HBOXES -1; # pos at right of -35 box **#_____*** 
			    $pos_tri_tmp[$sp][$i][2] = $pos_tri[$sp][$i][1];   # pos at left of -10 box ***_____#**
			    $pos_tri_tmp[$sp][$i][3] = $pos_tri[$sp][$i][1]+ $LOCAL_WIDTH_OF_R_HBOXES -1; # pos at right of -10 box ***_____**# 
			    $tab_lettre[$sp][$i] = ($H_l{'n'})x4;
			    # ($not_verbose)or print "tab_lettre $sp $i mis a 0\n";

			    # we get lower nd higher spacer between the two trinuc boxes
			    # my $sp_width = $pos_tri[$sp][$i][1] - $pos_tri[$sp][$i][0] - $WIDTH_OF_HBOXES;
			    # ($spacer_bounds[$sp][0] < $sp_width)or $spacer_bounds[$sp][0] = $sp_width;
			    # ($spacer_bounds[$sp][1] > $sp_width)or $spacer_bounds[$sp][1] = $sp_width;
			} 
			$sp > 1 and die "$prog_tag [Error] PB sp $sp >1 line ".__LINE__."\n";
		    }
		    if(! $not_verbose){

			print "verif pos_tri_tmp\n";
			foreach $sp(0..1)
			{
			    print "sp $sp\n";
			    foreach my $i(0..$#{$pos_tri_tmp[$sp]})
			    {
				print "pos_tri_tmp de $sp $i @{$pos_tri_tmp[$sp][$i]}\n";
			    }
			}
			print "fin verif pos_tri_tmp\n";
		     }


		    # OK die "sortedind 0 @{$sorted_ind[0]}\nsortedind 1 @{$sorted_ind[1]}\n";


		    # for  $sp(0..1){
		    for  ($sp = 0; $sp <= $#pos_tri; $sp++){

			($not_verbose)or print "\n**********************\nBacterie $sp\n";
			
			my $file_whole_genome;
			find sub { 
			  if($File::Find::name=~ /$deb_file_whole_genome$name_fic_sp[$sp]/){
			    $file_whole_genome = $File::Find::name; 
			  }
			}, "$way_to_whole_genome";
			  open(FWG, '<', $file_whole_genome)or die "$prog_tag [Error] Can not open file_whole_genome $file_whole_genome file of $way_to_whole_genome dir:$!, line ".__LINE__."\n"; 
				    
			  # USED for statistics in sub RAPPORT_CRITERIA
			  # not_fasta_whole_seq0_sco.txt
			  # ($pid = open(LSWHOLE,'ls '.$way_to_whole_genome.$deb_file_whole_genome.$name_fic_sp[$sp]."* |"))or die "Impossible to open ls command line ".__LINE__.": $!";
			  # (kill 0, $pid) or die "ls invocation failed:$!";



			# table used to record matrices associated with motifs
			@display = ();

			# spacer_bounds store spacer limits used to generate consensus motif with scoring matrix (updated after)
			# we intialize it in RECUP_LETTRE_SORT
			my @border_words  = (lc($subw1), lc($subw2));
			my @proba_ext     = (); # only declared because we need a parameter (useless for first call)
			@spacer_bounds = ();
			@interv_int    = (); # will record limits of the intervals of sequences which are interesting according to the middle of sequences scores. Only these intervals will be displayed
			# record seq ind sorted according to xor results****************************************
			@sorted_ind = ();
			# initialization
			for my $i(0..$#{ $seq_sp[$sp] }){ push @sorted_ind, $i; };		    
			
			if($border_words[0] =~ /n/){ @bool_first_call[0..1] = (0,0); }
			else                       { @bool_first_call[0..1] = (1,1); }
			if($border_words[1] =~ /n/){ @bool_first_call[2..3] = (0,0); }
			else                       { @bool_first_call[2..3] = (1,1); }

			# name of motif directory where are located all computed ratio for motifs for a given min_nb_of_waited_sigma_motifs
			$MOTIF_DIR_bact_min_nb_of_waited_sigma_motifs = $MOTIF_DIR.$name_fic_sp[$sp].'_'.$min_nb_of_waited_sigma_motifs.'/';

			# we get number of fixed letters in trinuc
			# by getting trinuc lengths...
			# my $nb_of_l_after_extend_to_evaluate = length($border_words[0]) + length($border_words[1]);

			# and withdrawing(soustraction) of number of 'n' character found (avoid to store each nb of letter foreach trinuc...if we would like to compute only once at the beginning)
			# while($border_words[0] =~ /n/g){ $nb_of_l_after_extend_to_evaluate--; }
			# while($border_words[1] =~ /n/g){ $nb_of_l_after_extend_to_evaluate--; }

			$nb_of_l_after_extend_to_evaluate = $nbl[ substr($str_trinuc_give_VAR, $num_subw1, 1) ] + $nbl[ substr($str_trinuc_give_VAR, $num_subw2, 1) ];
			# ($subw2 =~ /ggta/)and print "nb_of_l_after_extend_to_evaluate $nb_of_l_after_extend_to_evaluate pour $border_words[0] $border_words[1] correspondant respectievement en nb de lettres a ".$nbl[ substr($str_trinuc_give_VAR, $num_subw1, 1) ].' et '.$nbl[ substr($str_trinuc_give_VAR, $num_subw2, 1) ]."\n";
			@initial_border_words = ($border_words[0], $border_words[1]);

			$initial_total_boxes_length     = length($border_words[0]) + length($border_words[1]);
			$initial_total_nb_boxes_letters = 0;

			while($border_words[0] =~ /[$alphabet_trinuc]/g){ $initial_total_nb_boxes_letters++; }
			while($border_words[1] =~ /[$alphabet_trinuc]/g){ $initial_total_nb_boxes_letters++; }


			($not_verbose)or print "INITIAL border_words vaut @initial_border_words, initial_total_nb_boxes_letters $initial_total_nb_boxes_letters, initial_total_boxes_length $initial_total_boxes_length\n";

			# for debug ?? *****
			$time_long_or_not = time();
			my $duree = $time_long_or_not - $prev_time_long_or_not;
			$not_verbose or print "************************************************\ntreatment of $subw1, $subw2, bact $sp, at ".scalar(localtime(time())).", duree du precedent $duree\n";
			# if($duree > 1800){
			  #  die "last treatment last $duree seconds!!\n";
			# }
			# else{
			    # $prev_time_long_or_not = $time_long_or_not;
			# }
			# end for debug ?? *****


			my $string_let_pos_p_main      = ''; 
			my $bool_not_recorded_interval = 1;

			# var for displaying of results for one interval of one species and score computing
			if($sp){
			  $fic_sp = \*FALIGNSP2; 
			  $handle_ref_glob_sp0 = $handle_ref_glob10;
			  $handle_ref_glob_sp1 = $handle_ref_glob11;
			}
			else{ 
			  $fic_sp = \*FALIGNSP1; 
			  $handle_ref_glob_sp0 = $handle_ref_glob00;
			  $handle_ref_glob_sp1 = $handle_ref_glob01;
			}

			$length_same_begin_init = 18 + $LOCAL_WIDTH_OF_L_HBOXES + $LOCAL_WIDTH_OF_R_HBOXES;

			&RECUP_LETTRE_SORT_FIRST_CALL(0, $#{$pos_tri_tmp[$sp]}, @bool_first_call, $bool_not_recorded_interval, \@border_words, \@proba_ext, 1, $string_let_pos_p_main, 0, 0);

			$not_verbose or print "END treatment of $subw1, $subw2, bact $sp, at ".scalar(localtime(time())).", duree du precedent $duree\n************************************************\n";
			print {$fic_sp} "END treatment of $subw1, $subw2, bact $sp, at ".scalar(localtime(time())).", duree du precedent $duree\n************************************************\n";

				   $sp > 1 and die "$prog_tag [Error] PB sp $sp > 1 line ".__LINE__."\n";

				   close FWG; # close file_whole_genome
		    }	# end for $sp(0..$#pos_tri)	
		
		    $sp--;

		    # memory free
		    %{$H_nbw{$nbw}{$num_subw1}{$num_subw2}} = ();
		    %{$trinuc1_trinuc2_spacer{$num_subw1}{$num_subw2}} = ();
		    #VERIFICATION LIBERATION MEMOIRE ??
		    #($not_verbose)or print "On SUPPRIME H_nbw $subw1 $subw2 present $nbw et r_trinuc1_trinuc2_spacer correspondant\n";
		    
		} # end for my $num_subw2(
	    } # end for my $num_subw1(
    }

    # FIN TRAITEMENT DE LOTS DE TRINUCS ****************************
    
    if($bool_print_sp_align){

	close(FALIGNSP1);
	close(FALIGNSP2);
    }

    undef %H_nbw;
    undef %trinuc1_trinuc2_spacer;

}

# close($handle_ref_glob00) or die "Impossible to close FSORSEQ$_, qui correspond a $handle_ref_glob00:$!";
# close($handle_ref_glob01) or die "Impossible to close FSORID$_, qui correspond a $handle_ref_glob01:$!";
# close($handle_ref_glob10) or die "Impossible to close FSORSEQ$_, qui correspond a $handle_ref_glob10:$!";
# close($handle_ref_glob11) or die "Impossible to close FSORID$_, qui correspond a $handle_ref_glob11:$!";
# `free`;
#system("echo \"\" | mail -s \"prog finished\" $author_mail");
# exit; 


# for my $c(@t_mots_surrep){
# print "$c $mots_surrep{$c}\n";
# }


# *****SUB TO COMPUTE SCORE ACCORDING TO MATRIX MOTIF**********************************
sub COUNT_SURROUNDING_RATES(\@$$$\@\$$)
{
    my ($r_stat, $deb_l, $end_l, $reverse, $r_nb_l_added_matrice, $r_tab_count, $r_limit_extension_interne, $interne) = @_;
    my $rate = 0.7;
    my $subword = '';
    my @translate = ('a', 'c', 'g', 't');

    ($not_verbose)or print "SUB COUNT_SURROUNDING_RATES @_\n";

    if($reverse)
    {
	# for each position
      BOX1:for my $j(reverse $deb_l..$end_l)
      {
	  # for each letter
	  my $bool_found = 0;
	  for my $l(0..$#{$r_stat->[$j]})
	  {
	      if($r_stat->[$j][$l] > $rate)
	      {
		  $$r_nb_l_added_matrice++;
		  # we add the letter to interesting word which will compose a part of the motif to search
		  $subword = $translate[$l].$subword;

		  # ($not_verbose)or print "add of $translate[$l] at the beg, subword $subword \n";
		  $bool_found = 1;
		  if($interne){
		      $r_limit_extension_interne--;
		      ($r_limit_extension_interne)or last BOX1;
		  }
		  next BOX1;
	      }
	  }
	  ($bool_found)or last;
	  # print "on retourne subword $subword\n";
	  # return $subword;
      }
    }
    else
    {
	# for each position
      BOX2:for my $j($deb_l..$end_l)
      {
	  # for each letter
	  my $bool_found = 0;
	  for my $l(0..$#{$r_stat->[$j]})
	  {
	      if($r_stat->[$j][$l] > $rate)
	      {
		  $$r_nb_l_added_matrice++;
		  # we add the letter to interesting word which will compose a part of the motif to search
		  $subword = $subword.$translate[$l];
		  # ($initial_border_words[1] =~ /ggta/)and print "add of $translate[$l] at the beg, subword $subword, r_nb_l_added_matrice incremented $$r_nb_l_added_matrice, j $j, l $l line ".__LINE__."\n";
		  # ($not_verbose)or print "add of $translate[$l] at the end, subword $subword \n";
		  $bool_found = 1;
		  if($interne){
		      $r_limit_extension_interne--;
		      ($r_limit_extension_interne)or last BOX2;
		  }
		  next BOX2;
	      }
	  }
	  ($bool_found)or last;
      }
    }
    # ($not_verbose)or print "on retourne subword $subword et nb_l_added_matrice vaut $$r_nb_l_added_matrice dans COUNT_SURROUNDING_RATES\n";
    return $subword;

}
# *****END SUB TO COMPUTE SCORE ACCORDING TO MATRIX MOTIF******************************




# param 0 $reg_exp
# param 1 ref a @difference
# VAR GLOB: param 2 ref a (ref) $r_pos_mot
# param 3 position du mot trouve (en cours)
# param 4 ref a %pos_mot_pour_ce_mot, permet de calculer les differences
# VAR GLOB: param 5 %words_for_a_given_spacer
# VAR GLOB: param 6 ref a (ref) $r_ID_diff_pos
# param 7 ref a (ref) $r_mots_surrep (mots surrepresentes dans le fichier de diff de MM)
# param 8 ID du gene ou le motif est trouve
# GLOBAL VAR %diff_trouvees_pour_ce_motif table de hachage dont les cles sont toutes les differences trouvees pour le motif en cours dans le fichier analyse
# param 10 ref a (ref) $r_nb_motifs_distincts_par_diff_par_fic table de hachage comptant pour chaque difference, le nombre de motifs distincts (pour un fichier)
sub CALCUL_DIFF($\@$$)
{
    my ($regexp, $r_difference, $pos_mot_en_cours_sub, $gene_IDs_sub) = @_;
    
    my $IDseq_match    = '';
    my $indice_tab_diff_avant_jout = $#$r_difference;
    
    # ($not_verbose)or print "SUB CALCUL_DIFF @_\n";

    # pour chacun des identifiants de seq correspondant aux positions enregistrees pour ce mot...
    while( my ($clef) = each %pos_mot_pour_ce_mot)
    {
	# s'il ne correspond pas a l'identifiant de la sequence dans laquelle on a trouve le motif en cours...
	if($clef !~ /$gene_IDs_sub/)
	{
	    # pour chacune des positions enregistrees...
	    for my $differentes_pos(@{$pos_mot_pour_ce_mot{$clef}})
	    {
		# on calcule la difference de position entre la position du match trouve et la position de ce mot dans l'autre sequence
		push @$r_difference, $pos_mot_en_cours_sub - $differentes_pos;
		# on memorise l'identifiant de la sequence correspondante
		$IDseq_match .= $clef.' ';
		# on memorise la position du motif de la premiere sequence parcourue
		# FAUX	les positions sont utilisees pour les calculs de difference, il ne doivent concerner que l'ecart dans une seule sequence, PAS DANS LES DEUX!!!!	
		# (exists $pos_mot{$r_difference[$#$r_difference]}{$regexp})or push @{$pos_mot{$$r_difference[$#$r_difference]}{$regexp}}, $differentes_pos;
	    }
	}
    }
    
    # pour toutes les differences enregistrees
    for my $difference($indice_tab_diff_avant_jout + 1..$#$r_difference)
    {
	# if the same motif with the same difference already exists, we count it and we had the corresponding gene ID to the dedicated variable
	if(exists $words_for_a_given_spacer{$r_difference->[$difference]}{$regexp})
	{ 
	    $ID_diff_pos{$r_difference->[$difference]}{$regexp} .= $gene_IDs_sub.' '; 
	    $words_for_a_given_spacer{$r_difference->[$difference]}{$regexp}++;
	    # print "$regexp IDseq $ID_diff_pos{$$r_difference[$difference]}{$regexp}\ndifference $$r_difference[$difference]\n";
	}
	else
	{ 
	    # if it is the first time we find this motif, we record we found it twice (as there are 2 sequences implicated and we had the IDs of the sequences in the dedicated variable
	    $words_for_a_given_spacer{$r_difference->[$difference]}{$regexp} = 2;
	    $ID_diff_pos{$r_difference->[$difference]}{$regexp} = $IDseq_match.$gene_IDs_sub.' '; 
	    # print "$regexp IDseq $ID_diff_pos{$$r_difference[$difference]}{$regexp}\ndifference $$r_difference[$difference]\n";
	}
	# we record the motif position in current sequence
	push @{$pos_mot{$r_difference->[$difference]}{$regexp}}, $pos_mot_en_cours_sub;

	# comptage du nombre de mots trouves pour une difference de motif donnee
	if(! exists $diff_trouvees_pour_ce_motif{$r_difference->[$difference]})
	{
	    $diff_trouvees_pour_ce_motif{$r_difference->[$difference]} = 1;
	    
	    if(exists $nb_motifs_distincts_par_diff_par_fic{$r_difference->[$difference]})
	    {
		$nb_motifs_distincts_par_diff_par_fic{$r_difference->[$difference]}++; 
	    }
	    else
	    {
		$nb_motifs_distincts_par_diff_par_fic{$r_difference->[$difference]} = 1; 
	    }
	}
    }   
}

sub RECHERCHE($$\@)
{
    # mot recherche 

    # VAR GLOB: reference a la table de hachage des mots surrepresentes dans les seq intergen div par rapport a conv, %words_for_a_given_spacer (for &calcul_diff)
    # VAR GLOB: reference a la table de hachage referencant les mots pour une difference de position donnee (comparaison de seulement deux sequences), %ID_diff_pos (for &calcul_diff)
    # VAR GLOB: %pos_mot   table de hachage stockant les positions du mot pour une sequence donnee (for &calcul_diff)

    my ($reg_exp, $class_motif, $r_l_seq) = @_; 

    # ($not_verbose)or print "SUB RECHERCHE @_\n";
    
    my $bool_deja_repertorie_pour_cette_seq = 0;
    # longueur du mot recherche
    my $length_word   = length($reg_exp); 
    my $gene_IDs      = '';
    # sequence en cours (stockee au fur et a mesure de la lecture du fichier)
    my $sequence      = ''; 

    # GLOBAL VAR FOR CALCUL_DIFF
    # va stocker (pour un motif et un fichier, temporairement) toutes les differences de position entre deux sequences pour un mot donne
    my @difference                  = (); 
    my $pos_mot_en_cours            = 0;
    %pos_mot_pour_ce_mot            = ();    
    %diff_trouvees_pour_ce_motif   = (); # GLOBAL VAR (local globalisee en fait)

    seek FPROM, 0, 0;
    
    # tant qu'on est dans le fichier
    while (my $line = <FPROM>)
    {
		# treatment of first seq
		# si on commence a lire une autre sequence (son identifiant)
		if(( $line =~ /^>/)and($sequence ne ''))
		{
		    # tant qu'on trouve le mot
		    
		    # bool_deja_repertorie_pour_cette_seq :booleen permettant de ne compter qu'une occurence du motif par sequence
		    $bool_deja_repertorie_pour_cette_seq = 0;
	
		    while($sequence =~ m/$reg_exp/g)
		    {  
				# $pos_mot_en_cours = pos($sequence) - $length_word + 1;
				$pos_mot_en_cours = $-[0] + 1;
				
				# en cas de premiere occurrence, on incremente le compteur correspondant au nombre de sequences dans lesquelles on trouve le motif
		
				# OPTIMISER
				&CALCUL_DIFF($reg_exp,\@difference,$pos_mot_en_cours,$gene_IDs);
				
				($bool_deja_repertorie_pour_cette_seq)or $bool_deja_repertorie_pour_cette_seq = 1;
				push @{ $pos_mot_pour_ce_mot{$gene_IDs}}, $pos_mot_en_cours;
		    }
		    
		    # on memorise l'intitule de la sequence en cours (la deuxieme au moins)
		    $line =~ /:([^\s]+)/;
		    $gene_IDs = $1;
		    # print "gene_IDs $gene_IDs 1\n";
	
		    # we store seq length
		    ($#$r_l_seq == -1)and push @$r_l_seq, length($sequence);
		    # on efface la sequence precedente
		    $sequence = '';
		}
		elsif($line =~  /^>.+?:([^\s]+)/)
		{
		    $gene_IDs = $1;
		    # print "gene_IDs $gene_IDs 2\n";
		}
		elsif($line =~  /^>/)
		{
		    die "$prog_tag [Error] Problem regexp line $line not treated, line ".__LINE__."\n";
		}
		else{
		    $sequence .= $line;
		    chomp $sequence;
		}
    }

    # we store seq length
    ($#$r_l_seq == 0)and push @$r_l_seq, length($sequence);
    # $r_l_seq->[1] = length($sequence);


    # treatment of second seq
    $bool_deja_repertorie_pour_cette_seq = 0;
    while($sequence =~ m/$reg_exp/g)
    {  
		# $pos_mot_en_cours = pos($sequence) - $length_word + 1; 
		$pos_mot_en_cours = $-[0] + 1;
	
		&CALCUL_DIFF($reg_exp,\@difference,$pos_mot_en_cours,$gene_IDs);
	
		($bool_deja_repertorie_pour_cette_seq)or $bool_deja_repertorie_pour_cette_seq = 1;
		push @{ $pos_mot_pour_ce_mot{$gene_IDs}}, $pos_mot_en_cours;
    }
}

# UTILITIES

# ************************************************************************************************ 
# get part of sequences found in a fasta file
sub GET2SEQFASTAMINMAX(\@$\@$$$$$\@\%)
{
    my ( $r_H_important_motifs_num, $ID_fic_l, $r_specific_important_motifs, $b_min1, $longueur1, $b_min2, $longueur2 , $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS, $r_infos_div, $r_prev_infos_div ) = @_;

    my $bool_seq1_not_already_rec = 1;
    my $bool_seq2_not_already_rec = 1;

    # ($not_verbose)or print "SUB GET2SEQFASTAMINMAX @_\n";

    # print "GET2SEQFASTAMINMAX ".scalar(keys %$r_prev_infos_div)." tab @$r_infos_div[0..1]\n";
    # print "infos_div 0 $r_infos_div->[0] 1 $r_infos_div->[1]\n";
    # ($not_verbose)or print "infos_div 0 1 @$r_infos_div[0..1]\n";
    # my $pack_r_infos_div_0_1 = 
    if(exists $r_prev_infos_div->{ "@$r_infos_div[0..1]" } ){
	$bool_seq1_not_already_rec = 0;
    }

    # print "infos_div 2 $r_infos_div->[2] 3 $r_infos_div->[3]\n";
    # ($not_verbose)or print "infos_div 2 3 @$r_infos_div[2..3]\n";
    if(exists $r_prev_infos_div->{ "@$r_infos_div[2..3]" } ){
	$bool_seq2_not_already_rec = 0;
    }

    $r_prev_infos_div->{ "@$r_infos_div[0..1]" } = 1;
    $r_prev_infos_div->{ "@$r_infos_div[2..3]" } = 1;

    my $decal1 = 0; # decallage dans la position de la boite -35 du a sa position en debut de seq (ne permet pas d'avoir 10 nt a gauche du trinuc dans notre morceau de seq
    my $decal2 = 0; # comme decal1 mais pour la sequence de l'autre bacterie

    # b_min1 correspond a la position du debut du morceau de sequence que nous utilisons (10 nt en amont de la boite -35)
    # si b_min1 est negatif, cela veut dire que le trinuc de la boite -35 est en debut de sequence, il faut donc modifier la longueur du morceau de sequence (on retranche le nombre de nt manquant au debut
    if(($b_min1 < 0)and($bool_seq1_not_already_rec))
    {
	# on diminue la longueur du morceau de sequence du nombre de nucleotides (nt) manquant au debut de la sequence reelle (nous aurons donc moins de 10 nt pour l'extension a gauche du trinucleotide de la boite -35)
	$longueur1 += $b_min1; 
	# enregistre le decallage occasionne pour la position du trinucleotide de la boite -35 au sein de notre morceau de seq
	$decal1 = $b_min1; 

	# print "b_min valait $b_min1 qu'on ajoute à longueur $longueur1, bmin vaut maintenant 0\n";

	# on enregistre le fait que notre morceau de sequence correspond au tout debut de la sequence reelle
	$b_min1 = 0;
    }
    # comme b_min1 mais pour la sequence de l'autre bacterie
    if(($b_min2 < 0)and($bool_seq2_not_already_rec))
    {
	# on diminue la longueur du morceau de sequence du nombre de nucleotides (nt) manquant au debut de la sequence reelle (nous aurons donc moins de 10 nt pour l'extension a gauche du trinucleotide de la boite -35)
	$longueur2 += $b_min2; 
	# enregistre le decallage occasionne pour la position du trinucleotide de la boite -35 au sein de notre morceau de seq
	$decal2 = $b_min2; 

	# print "b_min valait $b_min2 qu'on ajoute à longueur $longueur1, bmin vaut maintenant 0\n";

	# on enregistre le fait que notre morceau de sequence correspond au tout debut de la sequence reelle
	$b_min2 = 0;
    }

    open(FILE,"< $tab_ID_fic[ unpack('n',$r_H_important_motifs_num->[0][$ID_fic_l]) ][0]")or die "$prog_tag [Error] Impossible to open $tab_ID_fic[ unpack('n',$r_H_important_motifs_num->[0][$ID_fic_l]) ][0] file: $!, line ".__LINE__."\n";
    # open(FILE,"< $tab_ID_fic[ unpack('n',$$r_file) ][0]")or die "Impossible to open $tab_ID_fic[ unpack('n',$$r_file) ][0] file: $!";

    my $sequence = '';
    my $diff_l1 = 0;
    my $bool_sec_entete_found = 0;

    my ($shiftsp0,$shiftsp1) = unpack('cc', $r_H_important_motifs_num->[3][$ID_fic_l]);

    ($not_verbose)or print "we get shiftsp: $shiftsp0,$shiftsp1\n";

    while( my $line = <FILE>)
    {
	if( $line =~ /^\s*$/ )
	{
	    next;
	}
	elsif( $line =~ /^>/ )
	{
	    $bool_sec_entete_found++;
	    if( $sequence ne '')
	    {
		($bool_seq1_not_already_rec)or next;

		my $box35_position = unpack('n', $r_H_important_motifs_num->[2][$ID_fic_l]) -1;

		# uppercase for words given by RMES (-35 and -10 boxes)
		# substr(sequence, position of -35 box, length of -35 box)
		substr($sequence, $box35_position, length($r_specific_important_motifs->[0])) = uc(substr($sequence, $box35_position, length($r_specific_important_motifs->[0])));
		# substr(sequence, position of -35 box + length of -35 box + length of spacer + spacer shift for this sequence (= position of -10 box), length of -10 box)
		substr($sequence, $box35_position + length($r_specific_important_motifs->[0])+ unpack('U',$r_specific_important_motifs->[1])+ $shiftsp0 - $shiftsp1, length($r_specific_important_motifs->[2])) = uc(substr($sequence, $box35_position + length($r_specific_important_motifs->[0])+ unpack('U',$r_specific_important_motifs->[1])+ $shiftsp0-$shiftsp1, length($r_specific_important_motifs->[2])));

		if(! $not_verbose)
		{
		    print "IDfic $ID_fic_l\n";
		    print "r_H_important_motifs_num 2 IDfic (box 35 pos): $box35_position\n";
		    print "r_specific_important_motifs 0: ".$r_specific_important_motifs->[0]."\n";
		    print "r_specific_important_motifs 1: ".unpack('U',$r_specific_important_motifs->[1])."\n";
		    print "r_specific_important_motifs 2: ".$r_specific_important_motifs->[2]."\n";
		    print "shiftsp0 $shiftsp0\n";
		    print "$sequence\n";
		}

		# check
		if(lc(substr($sequence, $box35_position + length($r_specific_important_motifs->[0])+ unpack('U',$r_specific_important_motifs->[1])+ $shiftsp0-$shiftsp1, length($r_specific_important_motifs->[2]))) ne $r_specific_important_motifs->[2]){
		    print substr($sequence, $box35_position + length($r_specific_important_motifs->[0])+ unpack('U',$r_specific_important_motifs->[1])+ $shiftsp0-$shiftsp1, length($r_specific_important_motifs->[2]))." ne ".$r_specific_important_motifs->[2]." line ".__LINE__." SUB GET!!!!!!!\n";
		    exit;
		     
		}

		# diff_l1 compte le nombre de nucleotides manquant a la sequence reelle pour avoir une marge de 10 nt a droite du trinucleotide de la boite -10 dans notre morceau de seq
		if(($diff_l1 = $b_min1 + $longueur1 - length($sequence)) > 0)
		{
		    # print ''.($b_min1 + $longueur1)." doit etre sup a ".(length($sequence))."\n";

		    # si le nombre de nt est insuffisant, on diminue la longueur de notre morceau de seq
		    $longueur1 -= $diff_l1;

		    # print "recalcul de longueur qui vaut $diff_l1\n";
		}
		else
		{
		    # sinon on met diff_l1 a zero puisqu'aucun decallage ne sera necessaire par la suite pour le calcul de la position du trinucleotide de la boite -10
		    $diff_l1 = 0;
		}
		push @{$seq_sp[0]}, substr($sequence, $b_min1, $longueur1);

		# print "ancienne version b0: on veut la partie qui commence a ".($longueur1 + $decal1 - $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS - $LOCAL_WIDTH_OF_L_HBOXES + $diff_l1)." soit longueur $longueur1 + decal $decal1 - localleft $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS - local width $LOCAL_WIDTH_OF_L_HBOXES + diff_l1 $diff_l1, de longueur $LOCAL_WIDTH_OF_L_HBOXES de la chaine\n$seq_sp[0][$#{$seq_sp[0]}] de longueur $longueur1 du fichier $ID_fic_l\nb_min vaut $b_min1\n";
		# print "nouvelle version b0: on veut la partie qui commence a $b_min1 de longueur $longueur1 du fichier $ID_fic_l\n";
		# print "on veut la partie qui commence a ".($longueur1 + $decal1 - $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS - $LOCAL_WIDTH_OF_HBOXES + $diff_l)." soit longueur $longueur1 + decal $decal1 - localleft $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS - local width $LOCAL_WIDTH_OF_HBOXES + diff_l $diff_l, de longueur $LOCAL_WIDTH_OF_HBOXES de la chaine\n$seq_sp[0][$#{$seq_sp[0]}] de longueur $longueur1 du fichier $r_H_important_motifs_num->[0][$ID_fic_l]\nb_min vaut $b_min1\n";
		# print "prem seq\n";
		my $pos_tri_35_1 = $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS + $decal1;
		my $pos_tri_10_1 = $longueur1 - $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS - $LOCAL_WIDTH_OF_R_HBOXES + $diff_l1;

		# we get associated IDs and positions
		push @{$ID_f[0]}, "@$r_infos_div[0..1] (@$r_infos_div[2..3])"; 
		# print "ID_f 0 vaut ${$ID_f[0]}[$#{$ID_f[0]}]\n";

		# MISE EN MAJUSCULES DES BOITES DE TRINUC (remplace par mise en maj des boites donnees par RMES
		# substr($seq_sp[0][$#{$seq_sp[0]}], $pos_tri_35_1, $LOCAL_WIDTH_OF_HBOXES) = uc(substr($seq_sp[0][$#{$seq_sp[0]}], $pos_tri_35_1, $LOCAL_WIDTH_OF_HBOXES));
		# substr($seq_sp[0][$#{$seq_sp[0]}], $pos_tri_10_1, $LOCAL_WIDTH_OF_HBOXES) = uc(substr($seq_sp[0][$#{$seq_sp[0]}], $pos_tri_10_1, $LOCAL_WIDTH_OF_HBOXES));

		push @{$pos_tri[0]}, [$pos_tri_35_1, $pos_tri_10_1];
		# print "pos @{${$pos_tri[0]}[$#{$pos_tri[0]}]}\n";
		# print  "${$seq_sp[0]}[ $#{$seq_sp[0]} ]\n";
		$sequence = '';

		($bool_seq2_not_already_rec)or last;
	    } # end if $sequence ne ''
	} # end if line =~ /^>/
	else
	{
	    if( ($bool_seq1_not_already_rec and ($bool_sec_entete_found == 1))or
		($bool_seq2_not_already_rec and ($bool_sec_entete_found == 2)) ){ chomp($sequence .= $line); }
	}
    } # end while my $line
    
    if(( $sequence ne '')and($bool_seq2_not_already_rec))
    {

	my $box35_position = unpack('n', $r_H_important_motifs_num->[1][$ID_fic_l]) -1;
	# uppercase for words given by RMES (-35 and -10 boxes)
	# substr(sequence, position of -35 box, length of -35 box)
	substr($sequence, $box35_position, length($r_specific_important_motifs->[0])) = uc(substr($sequence, $box35_position, length($r_specific_important_motifs->[0])));
	# substr(sequence, position of -35 box + length of -35 box + length of spacer + spacer shift for this sequence (= position of -10 box), length of -10 box)
	substr($sequence, $box35_position + length($r_specific_important_motifs->[0])+ unpack('U',$r_specific_important_motifs->[1]), length($r_specific_important_motifs->[2])) = uc(substr($sequence, $box35_position + length($r_specific_important_motifs->[0])+ unpack('U',$r_specific_important_motifs->[1]), length($r_specific_important_motifs->[2])));

	if(! $not_verbose)
	{
	    print "IDfic $ID_fic_l\n";
	    print "r_H_important_motifs_num 2 IDfic (box 35 pos): $box35_position\n";
	    print "r_specific_important_motifs 0: ".$r_specific_important_motifs->[0]."\n";
	    print "r_specific_important_motifs 1: ".unpack('U',$r_specific_important_motifs->[1])."\n";
	    print "r_specific_important_motifs 2: ".$r_specific_important_motifs->[2]."\n";
	    print "shiftsp1 $shiftsp1\n";
	    print "$sequence\n";
	}

	# check
	if(lc(substr($sequence, $box35_position + length($r_specific_important_motifs->[0])+ unpack('U',$r_specific_important_motifs->[1]), length($r_specific_important_motifs->[2]))) ne $r_specific_important_motifs->[2]){
	    print substr($sequence, $box35_position + length($r_specific_important_motifs->[0])+ unpack('U',$r_specific_important_motifs->[1]), length($r_specific_important_motifs->[2]))." ne ".$r_specific_important_motifs->[2]." line ".__LINE__." SUB GET!!!!!!!!\n";
	    exit; 
	}

	my $diff_l2 = 0;
	# diff_l2 compte le nombre de nucleotides manquant a la sequence reelle pour avoir une marge de 10 nt a droite du trinucleotide de la boite -10 dans notre morceau de seq
	if(($diff_l2 = $b_min2 + $longueur2 - length($sequence)) > 0)
	{
	    ($not_verbose)or print "b_min2 $b_min2 + longueur2 $longueur2 egale ".($b_min2 + $longueur2)." doit etre sup a ".(length($sequence))."\n";

	    # si le nombre de nt est insuffisant, on diminue la longueur de notre morceau de seq
	    $longueur2 -= $diff_l2;

	    ($not_verbose)or print "recalcul de longueur qui vaut $longueur2 apres suppression de $diff_l2\n";
	}	    
	else
	{
	    # print "diff_l2 $diff_l2 inf a 0\n";

	    # sinon on met diff_l2 a zero puisqu'aucun decallage ne sera necessaire par la suite pour le calcul de la position du trinucleotide de la boite -10
	    $diff_l2 = 0;
	}

	if(($longueur2 + $b_min2 > length($sequence))or($longueur2 < 0)or($b_min2 < 0))
	{
	    die "$prog_tag [Error] seq to short in $tab_ID_fic[ unpack('n',$r_H_important_motifs_num->[0][$ID_fic_l]) ][0] file\nWe want to get pos $b_min2 de longueur $longueur2, sequence is $sequence : \n$!, line ".__LINE__."\n";
	}

	push @{$seq_sp[1]}, substr($sequence, $b_min2, $longueur2);

	# print "ancienne version b1: on veut la partie qui commence a ".($LOCAL_LEFT_RIGHT_BOARDS_WIDTHS -1 + $decal2)." de longueur $LOCAL_WIDTH_OF_L_HBOXES de la chaine\n$seq_sp[1][$#{$seq_sp[1]}] de longueur $longueur2 du fichier $ID_fic_l\nb_min vaut $b_min2\n";
	# print "nouvelle version b1: on veut la partie qui commence a $b_min2 de longueur $longueur2 du fichier $ID_fic_l\n";
	# print "deuxieme seq\n";
	my $pos_tri_35_2 = $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS + $decal2;
	my $pos_tri_10_2 = $longueur2 - $LOCAL_LEFT_RIGHT_BOARDS_WIDTHS - $LOCAL_WIDTH_OF_R_HBOXES + $diff_l2;

	# we get associated IDs and positions
	push @{$ID_f[1]}, "@$r_infos_div[2..3] (@$r_infos_div[0..1])"; 
	# print "ID_f 1 vaut ${$ID_f[1]}[$#{$ID_f[1]}]\n";

	# substr($seq_sp[1][$#{$seq_sp[1]}], $pos_tri_35_2, $LOCAL_WIDTH_OF_HBOXES) = uc(substr($seq_sp[1][$#{$seq_sp[1]}], $pos_tri_35_2, $LOCAL_WIDTH_OF_HBOXES));
	# substr($seq_sp[1][$#{$seq_sp[1]}], $pos_tri_10_2, $LOCAL_WIDTH_OF_HBOXES) = uc(substr($seq_sp[1][$#{$seq_sp[1]}], $pos_tri_10_2, $LOCAL_WIDTH_OF_HBOXES));

	push @{$pos_tri[1]}, [$pos_tri_35_2, $pos_tri_10_2];
	# print "pos @{${$pos_tri[1]}[$#{$pos_tri[1]}]}\n";
	# print "${$seq_sp[1]}[ $#{$seq_sp[1]} ]\n";
    }
    close(FILE);
    # if($debug){ die; }
}

# ************************************************************************************************
# This program get letters surrounding trinucleotids boxes (according to given parameter, we can extend one or several letters
#  get one letters on the left and right of -35 and -10 trinucleotids according to given parameters
# param 0: r_seq_sp :ref to an array of part of all the sequences concerned by the current 
# trinucleotid dyad (the substring goes from 10 nucleotids before -35 trinucleotids to 10 nucleotids 
# after -10 trinucleotids max)
# param 1: r_sorted_ind :ref to a table storing the indice of the sequences of the previous table: 
# indices are sorted according to sequence similarities
# param 2: r_pos_tri_tmp : ref to a table storing for each seq the positions of letters flanking -35 
# and -10 trinucleotids and which we extend
# param 3: r_tab_lettre :ref to a table storing for each seq a binary code corresponding to flanking 
# letters of trinucleotid boxes
# param 4: r_H_l: ref to a hash table (keys = letter, value = binary code)
# param 5: r_score: ref to a table storing scores for each seq (here, only to be given to SORT_IND... 
# subprog)
# param 6 and 7: first and last indices of the sequences concerned by current extension (indices of
# r_sorted_ind table which give the real indices of the seq in r_seq_sp)
# bool_first_call : boolean to know if spacer intervals have been split (to avoid to have motifs like ggaat\w{9,18}gtt only because there is one seq with a spacer of 9)... 
sub RECUP_LETTRE_SORT($$$$$$\$\@\@\@$$$$$)
{
    my @ind_progress = @_[0..1];
    my @bool = @_[2..5];

    # bool0, 1, 2, 3 allow to know which position letters need to be extended for compairison
    my ( $bool_not_recorded_interval_l, $r_border_words, $r_proba_ext, $bool_split_in_record_not_done, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice) = @_[6..$#_];
   
    # print "sub RECUP_LETTRE_SORT indices @ind_progress at ",scalar(localtime(time()))," line ".__LINE__."\n";

    if(! $not_verbose)
    {
	print ' ' x length($string_let_pos_p);
	print "SUB RECUP_LETTRE_SORT @_\n";
    }

    # my $local_score = ;
    # si tous les booleens de choix d'extension sont a 0, nous savons que nous n'avons plus a etendre le motif
    if(! ($bool[0] or $bool[1] or $bool[2] or $bool[3]) )
    {
	($not_verbose)or print "we stop extension as no letter needs to be extended\n";
	return;
    } 
    
    if(!defined $ind_progress[0] ){ die "$prog_tag [Error] ind_progress de 0 non def $ind_progress[0], line ".__LINE__."\n";}
    if(!defined $ind_progress[1] ){ die "$prog_tag [Error] ind_progress de 1 non def $ind_progress[1], line ".__LINE__."\n";}
    # ($not_verbose)or print "RECUP_LETTRE_SORT APPELE @ind_progress, @bool, $bool_not_recorded_interval_l\n";

    for my $id_seq($ind_progress[0]..$ind_progress[1])
    {
	# print "idseq $id_seq\n";
	my $ind_of_tab_of_ind_to_modif = $sorted_ind[$id_seq];
	
	if(! defined $ind_of_tab_of_ind_to_modif){ die "$prog_tag [Error] non def line ".__LINE__."\n"; }
	# print "pos_tri_tmp $sp $ind_of_tab_of_ind_to_modif 0 vaut $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0]\n";
	if(!defined $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0])
	{
	    print "pos_tri_tmp de $ind_of_tab_of_ind_to_modif 0 non defini\n";
	    die "$prog_tag [Error] $!, line ".__LINE__."\n";
	}
	
	# CAUTION: we use elsif and not if because we extend only one letter at each step
	if($bool[0])
	{
	    # if we are not at the beginning of the sequence, we can decrease the indice of the letter which will be used
	    if($pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0] > 0)
	    {
		$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0]--; # on se place sur la position avant la boite -35
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 0, 1) = $H_l{ substr( $seq_sp[$sp][ $ind_of_tab_of_ind_to_modif ], $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0],1) };
	    }
	    # else, we have to indicate that no extension will be possible for this position
	    else
	    {
		# ($not_verbose)or print "bool0 mis a 0 dans RECUP_!!\n";
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 0, 1) = $H_l{ 'n' };
	    }
	}
	elsif($bool[1])
	{
	    # if the letter we want to extend don't have already been used by extension of neighbourged box
	    if($pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][1] < $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][2]-1 )
	    {
		$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][1]++; # on se place sur la position apres la boite -35
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 1, 1) = $H_l{ substr($seq_sp[$sp][ $ind_of_tab_of_ind_to_modif ],$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][1],1) };
	    }
	    else
	    {
		# ($not_verbose)or print "bool1 mis a 0 dans RECUP_!!\n";
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 1, 1) = $H_l{ 'n' };
	    }
	}
	elsif($bool[2])
	{
	    # if the letter we want to extend don't have already been used by extension of neighbourged box
	    if($pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][2] > $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][1]+1)
	    {
		$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][2]--; # on se place sur la position avant la boite -10
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 2, 1) = $H_l{ substr($seq_sp[$sp][ $ind_of_tab_of_ind_to_modif ],$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][2],1) };
	    }
	    else
	    {
		# ($not_verbose)or print "bool2 mis a 0 dans RECUP_!!\n";
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 2, 1) = $H_l{ 'n' };
	    }
	}
	else # if($bool[3])
	{
	    # if we are not at le last letter of the seq, we can extend 
	    if($pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][3] < length( $seq_sp[$sp][$ind_of_tab_of_ind_to_modif] ) - 2)
	    {
		$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][3]++; # on se place sur la position apres la boite -10
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 3, 1) = $H_l{ substr($seq_sp[$sp][ $ind_of_tab_of_ind_to_modif ],$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][3],1) };
	    }
	    else
	    {
		# ($not_verbose)or print "bool3 mis a 0 dans RECUP_!!\n";
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 3, 1) = $H_l{ 'n' };
	    }
	}
    }
    
     $not_verbose or print "****************\nverif r_tab_lettre, bool vaut @bool,  line ".__LINE__."\n";
     for my $id_seq($ind_progress[0]..$ind_progress[1])
     {
	 my $lettres = $rev_H_l[ substr($tab_lettre[$sp][ $sorted_ind[$id_seq] ], 0,1) ].' '.$rev_H_l[ substr($tab_lettre[$sp][ $sorted_ind[$id_seq] ], 1,1) ].' '.$rev_H_l[ substr($tab_lettre[$sp][ $sorted_ind[$id_seq] ], 2,1) ].' '.$rev_H_l[ substr($tab_lettre[$sp][ $sorted_ind[$id_seq] ], 3,1) ].' ';
        $not_verbose or print "ind $id_seq, realind $sorted_ind[$id_seq]: $lettres".unpack('B*',$tab_lettre[$sp][ $sorted_ind[$id_seq] ]).'  '.$seq_sp[$sp][ $sorted_ind[$id_seq] ]."\n";
     }

    if(! $not_verbose)
    {
      print "****************\n";
	print ' ' x length($string_let_pos_p);
	print "FIN SUB RECUP_LETTRE_SORT @_\n";
    }

    # ($not_verbose)or print "Dans RECUP r_border_words vaut @$r_border_words\n";
    if($ind_progress[1] - $ind_progress[0] > 0)
    {
	# &INIT_PROBA_TAB($r_proba_ext, $r_border_words, $ind_progress[1] - $ind_progress[0], @bool);
	($not_verbose)or print "lancement SORT $ind_progress[0], $ind_progress[1], @bool , $bool_not_recorded_interval_l, $r_border_words, $r_proba_ext, $bool_split_in_record_not_done, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice:line ".__LINE__."\n";

	&SORT_IND($ind_progress[0], 
		  $ind_progress[1], 
		  @bool , 
		  $bool_not_recorded_interval_l, 
		  $r_border_words, 
		  $r_proba_ext, 
		  $bool_split_in_record_not_done, 
		  $string_let_pos_p, 
		  $nb_of_intern_l_added_to_boxes, 
		  $nb_l_added_matrice);
    }
    
    
}


sub RECUP_LETTRE_SORT_FIRST_CALL($$$$$$\$\@\@\@$$$$)
{
    my @ind_progress = @_[0..1];

    $not_verbose or print "sub RECUP_LETTRE_SORT_FIRST_CALL ind @ind_progress at ",scalar(localtime(time()))," line ".__LINE__."\n";
    my @bool = @_[2..5];
    # bool0, 1, 2, 3 allow to know which position letters need to be extended for compairison
    my ( $bool_not_recorded_interval_l, $r_border_words, $r_proba_ext, $bool_split_in_record_not_done, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice) = @_[6..$#_];

    # ($not_verbose)or print "SUB RECUP_LETTRE_SORT @_\n";
    if(! $not_verbose)
    {
	print ' ' x length($string_let_pos_p);
	print "SUB RECUP_LETTRE_SORT @_\n";
    }
    # ($initial_border_words[1] =~ /ggta/)and print "SUB RECUP_LETTRE_SORT_FIRST_CALL @_, sp $sp\n";

    # my $local_score = ;
    # si tous les booleens de choix d'extension sont a 0, nous savons que nous n'avons plus a etendre le motif
    
    if(!defined $ind_progress[0] ){ die "$prog_tag [Error] ind_progress de 0 non def $ind_progress[0], line ".__LINE__."\n";}
    if(!defined $ind_progress[1] ){ die "$prog_tag [Error] ind_progress de 1 non def $ind_progress[1], line ".__LINE__."\n";}
    # ($not_verbose)or print "RECUP_LETTRE_SORT APPELE @ind_progress, @bool , $bool_not_recorded_interval_l\n";

    for my $id_seq($ind_progress[0]..$ind_progress[1])
    {
	# print "idseq $id_seq\n";
	my $ind_of_tab_of_ind_to_modif = $sorted_ind[$id_seq];
	
	if(! defined $ind_of_tab_of_ind_to_modif){ die "$prog_tag [Error] non def line ".__LINE__."\n"; }
	# print "pos_tri_tmp $sp $ind_of_tab_of_ind_to_modif 0 vaut $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0]\n";
	if(!defined $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0])
	{
	    print "pos_tri_tmp de $ind_of_tab_of_ind_to_modif 0 non defini line ".__LINE__."\n";
	    die "$prog_tag [Error] pos_tri_tmp de $ind_of_tab_of_ind_to_modif 0 non defini line ".__LINE__.": $!\n";
	}
	
	if($bool[0])
	{
	    # if we are not at the beginning of the sequence, we can decrease the indice of the letter which will be used
	    if($pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0] > 0)
	    {
		$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0]--; # on se place sur la position avant la boite -35
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 0, 1) = $H_l{ substr( $seq_sp[$sp][ $ind_of_tab_of_ind_to_modif ], $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][0],1) };
	    }
	    # else, we have to indicate that no extension will be possible for this position
	    else
	    {
		# ($not_verbose)or print "bool0 mis a 0 dans RECUP_!!\n";
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 0, 1) = $H_l{ 'n' };
	    }
	}
	
	if($bool[1])
	{
	    # if the letter we want to extend don't have already been used by extension of neighbourged box
	    if($pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][1] < $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][2]-1 )
	    {
		$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][1]++; # on se place sur la position apres la boite -35
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 1, 1) = $H_l{ substr($seq_sp[$sp][ $ind_of_tab_of_ind_to_modif ],$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][1],1) };
	    }
	    else
	    {
		# ($not_verbose)or print "bool1 mis a 0 dans RECUP_!!\n";
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 1, 1) = $H_l{ 'n' };
	    }
	}
	
	if($bool[2])
	{
	    # if the letter we want to extend don't have already been used by extension of neighbourged box
	    if($pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][2] > $pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][1]+1)
	    {
		$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][2]--; # on se place sur la position avant la boite -10
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 2, 1) = $H_l{ substr($seq_sp[$sp][ $ind_of_tab_of_ind_to_modif ],$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][2],1) };
	    }
	    else
	    {
		# ($not_verbose)or print "bool2 mis a 0 dans RECUP_!!\n";
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 2, 1) = $H_l{ 'n' };
	    }
	}
	
	if($bool[3])
	{
	    # if we are not at le last letter of the seq, we can extend 
	    if($pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][3] < length( $seq_sp[$sp][$ind_of_tab_of_ind_to_modif] ) - 2)
	    {
		$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][3]++; # on se place sur la position apres la boite -10
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 3, 1) = $H_l{ substr($seq_sp[$sp][ $ind_of_tab_of_ind_to_modif ],$pos_tri_tmp[$sp][ $ind_of_tab_of_ind_to_modif ][3],1) };
	    }
	    else
	    {
		# ($not_verbose)or print "bool3 mis a 0 dans RECUP_!!\n";
		substr($tab_lettre[$sp][ $ind_of_tab_of_ind_to_modif ], 3, 1) = $H_l{ 'n' };
	    }
	}
    }

    if(! $not_verbose)
    {
    
      print "****************\nverif r_tab_lettre, bool vaut @bool,  line ".__LINE__."\n";
      for my $id_seq($ind_progress[0]..$ind_progress[1])
	{
	  my $lettres = $rev_H_l[ substr($tab_lettre[$sp][ $sorted_ind[$id_seq] ], 0,1) ].' '.$rev_H_l[ substr($tab_lettre[$sp][ $sorted_ind[$id_seq] ], 1,1) ].' '.$rev_H_l[ substr($tab_lettre[$sp][ $sorted_ind[$id_seq] ], 2,1) ].' '.$rev_H_l[ substr($tab_lettre[$sp][ $sorted_ind[$id_seq] ], 3,1) ].' ';
	  print "ind $id_seq, realind $sorted_ind[$id_seq]: $lettres".unpack('B*',$tab_lettre[$sp][ $sorted_ind[$id_seq] ]).'  '.$seq_sp[$sp][ $sorted_ind[$id_seq] ]."\n";
	}
      print "****************\n";
      
      print ' ' x length($string_let_pos_p);
      print "FIN SUB RECUP_LETTRE_SORT @_\n";
    }

    # ($not_verbose)or print "Dans RECUP r_border_words vaut @$r_border_words\n";
    if($ind_progress[1] - $ind_progress[0] > 0)
    {
	($not_verbose)or print "On lance SORT avec indices  $ind_progress[0], $ind_progress[1] dans sub RECUP, bool_not_recorded_interval_l $bool_not_recorded_interval_l, bool_split_in_record_not_done $bool_split_in_record_not_done, string_let_pos_p $string_let_pos_p, nb_of_intern_l_added_to_boxes $nb_of_intern_l_added_to_boxes, nb_l_added_matrice $nb_l_added_matrice; line ".__LINE__."\n";
	# &INIT_PROBA_TAB($r_proba_ext, $r_border_words, $ind_progress[1] - $ind_progress[0], @bool);

	&SORT_IND($ind_progress[0], 
		  $ind_progress[1], 
		  @bool , 
		  $bool_not_recorded_interval_l, 
		  $r_border_words, 
		  $r_proba_ext, 
		  $bool_split_in_record_not_done, 
		  $string_let_pos_p, 
		  $nb_of_intern_l_added_to_boxes, 
		  $nb_l_added_matrice);
    }
    
    # verif sort ind
    # for my $ind(0..$#$r_sorted_ind)
    # {
    # print "ind $ind: $seq_sp[$sp][ $sorted_ind[$ind] ]\n";
    # }
    # die "affichage seq fini\n";
   
}


# ***********************************************************************************************

# GLOBAL VAR : sort_ind_to_sort : ref to a table of real seq indices (THIS table is sorted, not realy the seq)
# ind_beg_to_sort    : (relative) indice in sort_ind_to_sort of the first seq to sort
# ind_end_to_sort    : (relative) indice in sort_ind_to_sort of the last seq to sort
# GLOBAL VAR : score            : ref to score table (score for a seq corresponds to letters present in extension used)
# GLOBAL VAR : tab_lettre_sp    : ref to table of bordering letters foreach seq
# GLOBAL VAR : seq_sp         : ref to table of real seq (only for the concerned bacteria)
# GLOBAL VAR : pos_tri_tmp    : ref to table of positions of letters to extend foreach seq
# GLOBAL VAR : H_l            : ref to H table which gives numeric cod for corresponding letter
# bool                          : 4 booleans to know which position can be "extended"
# GLOBAL VAR : interv_int       : table storing first and last relative indices of sequences of an interesting set
# bool_not_recorded_interval_l  : boolean to know if a motif has already been recorded for this seq set
# GLOBAL VAR : spacer_bounds    : table corresponding to lower and higher spacer for the consensus motif corresponding to consensus motif
# GLOBAL VAR : display          : ref to table used for interesting matrices displaying
# GLOBAL VAR : name_fic_sp      : ref to table of species ID (for result files naming)
# GLOBAL VAR : way_to_whole_genome : ref to string corresponding to way_to_whole_genome files
# GLOBAL VAR : rapports         : ref to minimal threshold for R
# GLOBAL VAR : filtre_nb_l_added_matrice_maxi : ref to maximal threshold for extensio length
# GLOBAL VAR : pos_tri          : ref to table of INITIAL positions of letters to extend foreach seq (used to create matrice)
# r_border_words                : ref to table of trinucleotids used for extension proba
# r_proba_ext                   : ref to table of computed proba for extension (avoid to recompute at each step if the trinuc are unchanged)
# GLOBAL VAR : ID_f             : ref to table of files ID (so, genes ID)
# GLOBAL VAR : num_b            : number of the treated bacteria (0 for the first, 1 for the second)
# GLOBAL VAR : nb_upstr_seq     : number of upstream sequences considered :o)
# GLOBAL VAR : upstr_nt         : ref to nb of nucleotids in upstream seq
# GLOBAL VAR : total_nt         : ref to nb of nucleotids in total genome (plasmid also)
# GLOBAL VAR : spacer_shift_fic : table of authorized shifts in spacer between two files of upstream seq of orthologues
# bool_split_in_record_not_done : boolean used to know if spacer intervals have already been split (to avoid to do several times the same operation)
# string_let_pos_p                    : string to know letter en related proba for interesting motifs
# ************************************************************************************************
# This program sort seq in an interval according to letters at given positions

sub SORT_IND($$$$$$\$\@\@$$$){

    my ($ind_beg_to_sort, $ind_end_to_sort) = @_[0..1];
    my @bool = @_[2..5];
    my ( $bool_not_recorded_interval, $r_border_words, $r_proba_ext, $bool_split_in_record_not_done, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice) = @_[6..$#_];
    
    # avoid to erase for other call not related to same motif
    my @border_words = @$r_border_words; 

    # ($not_verbose)or print "SUB SORT_IND @_\n";

    if(! $not_verbose)
    {
      print "sub SORT at ",scalar(localtime(time())),", ind_beg_to_sort, ind_end_to_sort: $ind_beg_to_sort, $ind_end_to_sort\n";
	print ' ' x length($string_let_pos_p);
	print "SUB SORT_IND @_\n";
    }

    # ($not_verbose)or print "Appel SORT $r_sort_ind_to_sort, $ind_beg_to_sort, $ind_end_to_sort, $bool0, $bool1, $bool2, $bool3 , $bool_not_recorded_interval_l\n";
     
    ($ind_end_to_sort - $ind_beg_to_sort +1 < $nb_mini_of_seq_by_interv)and return;

    # count initialized to one because we have to count the initial seq other seq are compaired to
    my @global_cpt_score = (); 

    ($not_verbose)or print "ind_beg_to_sort $ind_beg_to_sort ind_end_to_sort $ind_end_to_sort\n";

    # allows to keep an interval only if DIFFERENT sequences are implicated (typically for case where one box is a polyN)
    my %ID_count    = ();
    my $count_ID    = 0;
    my @count_bad_letters_by_pos = ();

    for my $id_seq($ind_beg_to_sort..$ind_end_to_sort){

	my $real_ind_seq = $sorted_ind[$id_seq];
	
	# count of number of each letter in surrounding positions (indirect, by storing relatives positions)
	for my $p(0..3){

	    # $l == 4 means no letter, end of strin reached
	    # lowercase required to consider lower and uppercase equals (a eq A), and the letter can be used for trinuc composition 
	    my $l = substr($tab_lettre[$sp][ $real_ind_seq ], $p, 1);
	   # if(($l == 0)or($l == 3)){ print "Pour LOW_PROBA plus tard: p vaut $p on met $id_seq dans global_cpt_score $p ".substr($tab_lettre[$sp][ $real_ind_seq ], $p, 1).", seq $seq_sp[$sp][$real_ind_seq], $real_ind_seq, id_seq $id_seq, sp $sp\n"; }
	    if($l == 4){ 
	      push @{ $count_bad_letters_by_pos[$p] }, $real_ind_seq;    
	    }
	    else{ 
	      push @{ $global_cpt_score[$p][ $l ] }, $real_ind_seq; 
	    }
	    # print "Pour LOW_PROBA plus tard: p vaut $p on met $id_seq dans global_cpt_score $p ".substr($tab_lettre[$sp][ $real_ind_seq ], $p, 1).", seq $seq_sp[$sp][$real_ind_seq], $real_ind_seq, id_seq $id_seq, sp $sp\n"; 
	}
	if((! exists $ID_count{ $real_ind_seq })and($count_ID < $nb_mini_occ))
	{
	    # ($not_verbose)or print "nouvelle seq distincte $1\n";
	    $ID_count{ $real_ind_seq } = 1;
	    $count_ID++;
	}

    }
    # ($not_verbose)or print "\n";

    # verif count_bad_letters_by_pos OK
    # print "verif count_bad_letters_by_pos\n";
    # for my $posv(0..$#count_bad_letters_by_pos){
      # print "posv $posv, id:\n";
      # foreach(@{ $count_bad_letters_by_pos[$posv] }){
	# print "$_ ";
      # }
      # print "\n";
    # }

    # verif seq OK
    # print "SEQ IMPLIQUEES\n";
    # for(@$r_seq_sp_l)
    # {
    # print "$_\n";
    # }
    # print "\n";

    if(! $not_verbose){

      print "verif comptage OK\n";
      for my $i(0..3){
	
	print "pos $i: ";
	
	for my $l(0..$#{ $global_cpt_score[$i] }){

	  print "lettre $l ($rev_H_l[$l]) compte ".($#{ $global_cpt_score[$i][$l] } +1)." :@{ $global_cpt_score[$i][$l] }, ";
	}
	print "\n";
      }
      print "verif comptage END\n";
    }

    # on veut que le nombre de seq pour un intervalle interessant soit sup ou egale au nombre de seq mini par intervalle (= ensemble de seq)
    if($count_ID < $nb_mini_of_seq_by_interv )
    {
	if(! $not_verbose)
	{
	  print "nb_mini_of_seq_by_interv $nb_mini_of_seq_by_interv > nbseq interv $count_ID ou < $nb_mini_occ et bool_not_recorded_interval $bool_not_recorded_interval, line ".__LINE__."\n";
	    print ' ' x length($string_let_pos_p);
	    print "FIN SUB SORT_IND @_\n";
	}
	return;
    } # avoid to divid by zero for score computing (which is useless in this case)
    else{
	($not_verbose)or print "on continue car count_ID $count_ID >= $nb_mini_of_seq_by_interv, line ".__LINE__."\n";
    }

    &INIT_PROBA_TAB($r_proba_ext, 
		    \@border_words, 
		    $ind_end_to_sort-$ind_beg_to_sort);

    # if(! $not_verbose)
    # {
    # &verif_low_proba(\@low_proba, $pos_proba_min, __LINE__);
    # }


    # proba have to be computed to know which letter to extend 
    my $min = 1;
    my $nb_seq_max = 0; # max number of seq involved in quite interesting letter proba
    my $nb_seq_tmp = 0; # number of seq involved in min proba compute
    my $nb_seq     = 0; # number of seq involved with minimal proba
    # will be used to filter with affectation foreach interesting proba....

    my $pos_proba_min = -1; # store pos of lower proba if we take into account every position and every letter
    my @low_proba  = ();

    # just for information!!!! ??
    my $nb_total_of_low_proba = 0;

    for my $pos(0..3)
    {
	for my $l(0..$#{ $global_cpt_score[$pos] })
	{
	    # proba compute
	    my $proba = 0;

	    my $ind_l;
	    # if( $#{ $global_cpt_score[$pos][$l] }+1 < $nb_mini_of_seq_by_interv )
	    if( $#{ $global_cpt_score[$pos][$l] }+1 < $nb_mini_occ )
	    {
		$proba = 1;
		# ($not_verbose)or print "$l letter is only present ".($#{$global_cpt_score[$pos][$l]}+1)." times, not suffisant (< $nb_mini_of_seq_by_interv)\n";
	    }
	    # if extension on a letter at left, $pos & 1 == 0
	    # if extension on a letter at right, $pos & 1 == 1 (proba different if left or right)
	    else
	    {
		# ${ $global_cpt_score[$pos]{$k} }[0]] give relative indice for seq
		# $sorted_ind[ ... ] give real indice for seq
		# $tab_lettre[$sp][ ... ] give the letters bordering trinuc
		# substr is used to extract interesting letter (coded by bits)
		$ind_l = substr( $tab_lettre[$sp][ $global_cpt_score[$pos][$l][0] ], $pos, 1);
		# ($not_verbose)or print "ind_l $ind_l ($rev_H_l[$ind_l]), pos vaut $pos \n";
		$nb_seq_tmp = $#{ $global_cpt_score[$pos][$l] }+1;
		# ($not_verbose)or print "r_border_words vaut @$r_border_words\n";
		foreach(@{ $r_proba_ext->[ $pos ][ $ind_l ] }[ $nb_seq_tmp..$#{ $r_proba_ext->[ $pos ][ $ind_l ] } ]){  $proba += $_; }
		# ($not_verbose)or print "proba $proba pour que la lettre $ind_l apparaisse au moins $nb_seq_tmp fois apres nos $#{ $r_proba_ext->[ $pos ][ $ind_l ] } trinuc\n";
	    }
	    
	    
	    # if( $#{ $global_cpt_score[$pos]{$k} }+1 > $max )
	    
	    if($proba <= $max_valid_prob_ext)
	    {
		push @{$low_proba[$pos]}, [$ind_l, $proba];
		# print "dans low_proba: lettre $rev_H_l[$ind_l], p: $proba, pos $pos, concerne $nb_seq_tmp seq\n";

		$nb_total_of_low_proba++;

		if( $proba  < $min )
		{
		    $nb_seq = $nb_seq_tmp;
		    $min = $proba;
		    
		    ($not_verbose)or print "MIN $min vaut $proba concerne $nb_seq_tmp seq\n";
		    $pos_proba_min = $pos;
		}
	    }

	    ($proba <= $max_valid_prob_ext)and($nb_seq_max < $nb_seq_tmp )and $nb_seq_max = $nb_seq_tmp;
	}
    }

    if(($pos_proba_min == -1)and($#low_proba != -1)){
      print "PB: pos_proba_min $pos_proba_min BUT low_proba tab last ind $#low_proba, nb_total_of_low_proba $nb_total_of_low_proba, min $min\n";
      for(0..3){ &verif_low_proba(\@low_proba, $_, __LINE__); }
      print "PB line ".__LINE__."\n";
      die "$prog_tag [Error] PB , see debug file tail!!, line ".__LINE__."\n";
    }
    if($pos_proba_min == -1){
      $not_verbose or print "RETURN car pos_proba_min $pos_proba_min\n";
      @low_proba = ();
      return 0;
    }

    # sort of every proba <= $max_valid_prob_ext
    # we sort not only position with lower proba because we are going to treat every low proba (with stack)

    # on place la lettre de plus faible proba au debut (on l'echange avec la place qu'elle occupe) 

    for(my $pos_i = 0; $pos_i < 4; $pos_i++){
      
      ($#{ $low_proba[$pos_i] } == -1)and next;

      $not_verbose or &verif_low_prob(\@low_proba, $pos_i, __LINE__);      

      (defined $low_proba[$pos_i])and @{$low_proba[$pos_i]} = sort { $a->[1] <=> $b->[1] } @{$low_proba[$pos_i]};
      
      $not_verbose or &verif_low_prob(\@low_proba, $pos_i, __LINE__);
    }
    
    
    my ($beg_loop, $end_loop, @const_mem_sorted_ind, $mem_bool_not_recorded_interval, @const_tab_lettre, @const_pos_tri_tmp); 
    if($bool_use_every_low_proba){ 
	($beg_loop,$end_loop) = (0,3);
	@const_mem_sorted_ind = @sorted_ind; 
	$mem_bool_not_recorded_interval = $bool_not_recorded_interval;
	# print "ori: tab_lettre @tab_lettre\n";
	@const_tab_lettre = @{ &clone(\@{ $tab_lettre[$sp] }) };
	# print "clone: @const_tab_lettre @tab_lettre\n";
	@const_pos_tri_tmp    = @{ &clone( $pos_tri_tmp[$sp] ) };
	$not_verbose or print "clone pos_tri_tmp line ".__LINE__."\n";
	($#{ $pos_tri_tmp[$sp] } != $#const_pos_tri_tmp)and die "$prog_tag [Error] clone does not work line ".__LINE__."\n";
    }
    else{ ($beg_loop,$end_loop) = ($pos_proba_min,$pos_proba_min); }


    # EN COURS ***********************************************************************************************

    # LOOP TO TREAT EVERY LOW PROBA, if beg_loop == end_loop, we treat every proba at pos proba_min, NOT EVERY **************************
    for(my $pos_i = $beg_loop; $pos_i <= $end_loop; $pos_i++){ 

	if((not defined $low_proba[$pos_i])or($#{ $low_proba[$pos_i] } == -1)){ next; }

	if($bool_use_every_low_proba and($pos_i != $beg_loop)){ 
	  # REINITIALIZATION of data TO TREAT OTHER POS of low_proba than proba_min
	    @sorted_ind = @const_mem_sorted_ind; 
	    @{ $tab_lettre[$sp] } = @{ &clone(\@const_tab_lettre) };
	    $pos_tri_tmp[$sp] = &clone(\@const_pos_tri_tmp );
	    $not_verbose or print "clone pos_tri_tmp line ".__LINE__."\n";
	    $bool_not_recorded_interval = $mem_bool_not_recorded_interval;
	    ($#{ $pos_tri_tmp[$sp] } != $#const_pos_tri_tmp)and die "$prog_tag [Error] clone does not work line ".__LINE__."\n";
	}

	# indice used to browse @{ $global_cpt_score[ $pos_i[0] ]{ $low_proba[ $pos_i ][ 0 ] } } table
	my $first_ind_replacement_tab   = 0; 
	
	# recuperation des limites de taille de spacer en profitant de la boucle ci-dessous de parcours des seq de l'intervalle
	my @tempo_spacer_bounds = @tempo_spacer_bounds_init; # globaly for all seq related to a low proba
	my @tempo_spacer_bounds_interv; # for each seq set corresponding to the pos with the lower proba
	
	# allows to keep an interval only if DIFFERENT sequences are implicated (typically for case where one box is a polyN)
	%ID_count       = ();
	my @count_ID    = ();
	
	my $exchanged_real_ind = '';

	# store last indice corresponding to the last good seq of an interesting seq set    
	my @last_ind_seq_set = ($ind_beg_to_sort-1); 
	# print "initialisation de last_ind_seq_set, dernier indice $#last_ind_seq_set, de valeur $last_ind_seq_set[$#last_ind_seq_set]\n";
	# push @last_ind_seq_set, -1; 
	
	# my $ind_int          = $ind_beg_to_sort; # for my $ind_int($ind_beg_to_sort..$ind_end_to_sort)
	# print "ind_int $ind_int initialise\n";
    
	# (defined $low_proba[$pos_i])or die "low_proba de pos_proba_min $pos_i not def line ".__LINE__."\n"; 
	# print "nb_lettre concerned by pos $pos_i: $#{ $low_proba[$pos_i] }+1\n";
	
	# print "seq AVANT TRI @sorted_ind[$ind_beg_to_sort..$ind_end_to_sort]\n";
	
	# foreach(@sorted_ind){ print "seq $seq_sp[$sp][$_]\n"; }
	
	# @l_sorted_by_proba contient les lettres par ordre croissant de proba
	
	
	my $relative_interest_ind = -1;
	my @new_sort_ind = ();
	my %treated_l = ();
	
	# TREATMENT OF 1 proba at THIS pos *********************************************************************
	for(my $ind_l_this_pos = 0; $ind_l_this_pos <= $#{ $low_proba[$pos_i] }; $ind_l_this_pos++){
	    
	    push @tempo_spacer_bounds_interv, [ @tempo_spacer_bounds_init ];
	    push @count_ID, 0;
	    %ID_count = ();
	    
	    
	    # we sort ind in new_sort to have seq corresponding to low proba in firsts position	    
	    push @new_sort_ind, @{ $global_cpt_score[$pos_i][ $low_proba[$pos_i][$ind_l_this_pos][0] ] };

	    # print "treatment of letters $low_proba[$pos_i][$ind_l_this_pos][0], concerns $#{ $global_cpt_score[$pos_i][ $low_proba[$pos_i][$ind_l_this_pos][0] ] } seq for pos $pos_i line ".__LINE__."\n";
	    # print "longueurs espaces\n";
	    $treated_l{ $low_proba[$pos_i][$ind_l_this_pos][0] } = 1;

	    for($first_ind_replacement_tab = 0;  $first_ind_replacement_tab <= $#{ $global_cpt_score[$pos_i][ $low_proba[$pos_i][$ind_l_this_pos][0] ] }; $first_ind_replacement_tab++){
		
	      # print "first_ind_replacement_tab $first_ind_replacement_tab correspond a low_proba @{ $global_cpt_score[$pos_i][ $low_proba[$pos_i][$ind_l_this_pos][0] ] } ln ".__LINE__."\n";
		
	      # we get real ind of the seq
	      $exchanged_real_ind = $global_cpt_score[$pos_i][ $low_proba[$pos_i][$ind_l_this_pos][0] ][ $first_ind_replacement_tab ]; # added 070306
		
	      
	      # count spacer intervalle
	      
	      # $pos_tri[$sp][$exchanged_real_ind][0]; # pos at right of -35 box #**_____*** 
	      # $pos_tri[$sp][$exchanged_real_ind][1]; # pos at left of -10 box ***_____#**
	      
	      
	      # we get length of spacer between trinuc for first seq
	      my $length = $pos_tri[$sp][$exchanged_real_ind][1] - $pos_tri[$sp][$exchanged_real_ind][0] - $LOCAL_WIDTH_OF_L_HBOXES;
	      # print "longueur espace $length = $pos_tri[$sp][$exchanged_real_ind][1] - $pos_tri[$sp][$exchanged_real_ind][0] - $LOCAL_WIDTH_OF_L_HBOXES\n"; # pour $seq[$sp][$sorted_ind]\n";
	      
	      if($length > 0){
		
		($tempo_spacer_bounds[0] < $length)or $tempo_spacer_bounds[0] = $length;
		($tempo_spacer_bounds[1] > $length)or $tempo_spacer_bounds[1] = $length;
		
		($tempo_spacer_bounds_interv[$#tempo_spacer_bounds_interv][0] < $length)or $tempo_spacer_bounds_interv[$#tempo_spacer_bounds_interv][0] = $length;
		($tempo_spacer_bounds_interv[$#tempo_spacer_bounds_interv][1] > $length)or $tempo_spacer_bounds_interv[$#tempo_spacer_bounds_interv][1] = $length;
	      }
	      else{
		
		($not_verbose)or print "length vaut $length = $pos_tri[$sp][$exchanged_real_ind][1] - $pos_tri[$sp][$exchanged_real_ind][0] - $LOCAL_WIDTH_OF_L_HBOXES\n";
		# trinuc are neirbourghs, no space, minimal spacer has to be set to 0 (default value of $length)
		$length = 0;
		$tempo_spacer_bounds[0] = $tempo_spacer_bounds_interv[$#tempo_spacer_bounds_interv][0] = $length;
		# ($tempo_spacer_bounds_interv[$#tempo_spacer_bounds_interv][0] < $length)or $tempo_spacer_bounds_interv[$#tempo_spacer_bounds_interv][0] = $length;
	      }
		
	      # COUNT OF NUMBER OF DIFFERENT SEQ****************************************************
	      
	      $ID_f[$sp][$exchanged_real_ind] =~ /^([^\s]+)/;
	      if((! exists $ID_count{ $1 })and($count_ID[$#count_ID] < $nb_mini_occ)){
		
		# ($not_verbose)or print "nouvelle seq distincte $1\n";
		$ID_count{ $1 } = 1;
		$count_ID[$#count_ID]++;
	      }
		
	      # COUNT OF NUMBER OF DIFFERENT SEQ****************************************************
	      
	    }
	    # print "\n";

	    # we record in last_ind_seq_set last indices in final sort which interest us:
	    # we use + $ind_beg_to_sort to avoid to sort part of table unchanged
	    push @last_ind_seq_set, $#new_sort_ind + $ind_beg_to_sort;
	    # push @last_ind_seq_set, $#new_sort_ind; # MODIFIED!!!!!
	    # print "ajout ds last_ind_seq_set, dernier indice $#last_ind_seq_set, de valeur $last_ind_seq_set[$#last_ind_seq_set]\n";
	}
	# END TREATMENT OF 1 proba at THIS pos *********************************************************************
	
	# to know until which seq indice seq have been used to generate matrice
	# $last_ind_seq_set[$#last_ind_seq_set] CORRESPONDS TO LAST INSTERESTING SEQ IND

	# my $der_interesting_ind_in_sort = $#new_sort_ind;

	# if($last_ind_seq_set[$#last_ind_seq_set] != $der_interesting_ind_in_sort){
	  # print "CAUTION: $der_interesting_ind_in_sort (der_interesting_ind_in_sort) DIFF de $last_ind_seq_set[$#last_ind_seq_set] (last_ind_seq_set de $#last_ind_seq_set)\n";
	# }

	# ADD in new_sort ind of interesting seq (linked to low_proba) (so present in %treated)
	for my $not_treated_l(0..3){
	    if((exists $treated_l{$not_treated_l})or(not defined $global_cpt_score[$pos_i][ $not_treated_l ])or($#{ $global_cpt_score[$pos_i][ $not_treated_l ]} == -1)){ next; }

	    # print "new sort last ind was $#new_sort_ind, and becomes ";
	    push @new_sort_ind, @{ $global_cpt_score[$pos_i][ $not_treated_l ] };
	    # print "$#new_sort_ind line ".__LINE__." because global_cpt_score has ".@{ $global_cpt_score[$pos_i][ $not_treated_l ] }." boxes\n";
	}
	# %treated_l = ();

	# ADD in new_sort seq ind wich have NO letter at this position (to short)
	if(defined $count_bad_letters_by_pos[$pos_i]){
	  # print "new sort last ind was $#new_sort_ind, and becomes ";
	  push @new_sort_ind, @{ $count_bad_letters_by_pos[$pos_i] };
	  # print "$#new_sort_ind line ".__LINE__." because count_bad_letters_by_pos has ".@{ $count_bad_letters_by_pos[$pos_i] }." boxes\n";
	}

	# CHECK correct run
	if($ind_end_to_sort - $ind_beg_to_sort != $#new_sort_ind ){ 
	  # @sorted_ind[$ind_end_to_sort..$ind_beg_to_sort] = sort @sorted_ind[$ind_end_to_sort..$ind_beg_to_sort];
	  @new_sort_ind = sort @new_sort_ind;
	  print "sorted_ind "; foreach($ind_end_to_sort..$ind_beg_to_sort){print "$sorted_ind[$_] "; } print "\nnew_sorted_ind @new_sort_ind\n";
	  die "$prog_tag [Error] $ind_end_to_sort - $ind_beg_to_sort (".($ind_end_to_sort - $ind_beg_to_sort).") != $#new_sort_ind line ".__LINE__."\n"; 
	}
	
	if(not $not_verbose){
	  if(($border_words[0] eq 'aacc')and($border_words[1] eq 'ccg')){
	    print "sorted ind AVANT\n";
	    foreach(@sorted_ind[$ind_beg_to_sort..$ind_end_to_sort]){ print "seq $seq_sp[$sp][$_], real ind $_\n"; }
	  }
	}

	# print "new_sort_ind\n@new_sort_ind\n";
	# if(@sorted_ind[$ind_beg_to_sort..$ind_end_to_sort] != @new_sort_ind){ die "size of sorted_ind de $ind_beg_to_sort a $ind_end_to_sort ".@sorted_ind[$ind_beg_to_sort..$ind_end_to_sort]." diff of size of new_sort_ind ".@new_sort_ind." line ".__LINE__."\n"; }

	# we record sort
	@sorted_ind[$ind_beg_to_sort..$ind_end_to_sort] = @new_sort_ind;
	
	if(not $not_verbose){
	  if(($border_words[0] eq 'aacc')and($border_words[1] eq 'ccg')){
	    print "sorted ind APRES\n";
	    foreach(@sorted_ind[$ind_beg_to_sort..$ind_end_to_sort]){ print "seq $seq_sp[$sp][$_], real ind $_\n"; }
	  }
	}

	%treated_l = ();
	# $der_interesting_ind_in_sort = $ind_beg_to_sort + $der_interesting_ind_in_sort;
	
	($not_verbose)or print "seq APRES TRI @sorted_ind[$ind_beg_to_sort..$ind_end_to_sort]\n";

	# foreach(@sorted_ind[$ind_beg_to_sort..$ind_end_to_sort]){ print "seq $seq_sp[$sp][$_]\n"; if(not defined $seq_sp[$sp][$_]){ die "seq_sp de $sp $_ non def line ".__LINE__."\n"; } }
	# print "\n";
	my $bool_not_recorded_interval_l = $bool_not_recorded_interval;
	
	# JUSQUE LA TRI OK
	
	# EVALUATION of EVERY seq with low proba at this position (so border letter of regexp can be "multiletter")
	if(($bool_not_recorded_interval_l)and($last_ind_seq_set[$#last_ind_seq_set] - $last_ind_seq_set[0] >= $nb_mini_occ)){
	    my $regexp_motif           = '';
	    my $nb_l_added_matrice_l = $nb_l_added_matrice;
	    # print "min $min, pos_proba_min $pos_i, last_ind_seq_set @last_ind_seq_set\n";
	    # we do not need because we have tested before if pos_proba_min == -1 160106
	    # if($#last_ind_seq_set > 0){
	    $not_verbose or print "MATRICE_CREATE2 ran line ".__LINE__." nb_seq ".($last_ind_seq_set[$#last_ind_seq_set] - $last_ind_seq_set[0])."\n";
	    &MATRICE_CREATE2(\@border_words, $low_proba[$pos_i], \$regexp_motif, \@tempo_spacer_bounds, $global_cpt_score[$pos_i], $last_ind_seq_set[$#last_ind_seq_set] - $last_ind_seq_set[0], $nb_of_intern_l_added_to_boxes, \$nb_l_added_matrice_l, $pos_i);
	    
	    # if we do not find any interesting thing here, we will not find anything
	    if($nb_l_added_matrice_l == $nb_l_added_matrice){
	      $not_verbose or print "No new letter in MATRICE_CREATE2, we return line ".__LINE__."!!!\n"; 
	      next;
	    }
	    print "regexp_motif $regexp_motif\n";

	    if(not $not_verbose){
	      print "On lance EVAL line ".__LINE__."\n";
	      
	      if($regexp_motif eq "[at]aacc\\w{17,19}ccg"){
		print "Seq involved:\n";
		for($last_ind_seq_set[0]..$last_ind_seq_set[$#last_ind_seq_set]){
		  print "seq $seq_sp[$sp][$_]\n";
		}
		print "Seq involved END\n";
		&verif_low_prob(\@low_proba, $pos_i, __LINE__);
	      }
	    }

	    &EVALUATION($ind_beg_to_sort, $last_ind_seq_set[$#last_ind_seq_set], 1, \@tempo_spacer_bounds, $pos_i, $r_proba_ext, \@border_words, $string_let_pos_p, $bool_not_recorded_interval_l, \$regexp_motif, \@low_proba, $global_cpt_score[$pos_i], $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice_l, \@last_ind_seq_set);
	}
	else{
	    $not_verbose or print "No stop bool_not_recorded_interval_l $bool_not_recorded_interval_l == 0, or $last_ind_seq_set[$#last_ind_seq_set] - $last_ind_seq_set[0] < $nb_mini_occ\n";
	    
	    # go to next position of low_proba if we treat every low_proba
	    next;
	}
	# if interesting motif, we look at 
	
	
	# print "nb_total_of_low_proba $nb_total_of_low_proba concerns ".($#low_proba+1)." positions\n";
	
	# FIN EN COURS ***********************************************************************************************
	
	
	
	# treatment of other seq
	if($#low_proba > 0){
	    # my @other_ind = ($last_ind_seq_set[$#last_ind_seq_set]+1, $ind_end_to_sort); 
	    my @other_ind = ($last_ind_seq_set[$#last_ind_seq_set]+1, $ind_end_to_sort); 
	    
	    if($other_ind[1] - $other_ind[0] +1  >= $nb_mini_occ)
	    {
		my $bool_local_not_recorded_interval = 0;
		# print "$other_ind[1] - $other_ind[0] +1  >= $nb_mini_of_seq_by_interv\n";	    
		# ($other_ind[1] > $r_interv_int->[$#$r_interv_ind][1]) means sort concerned sequences involved in an already recorded interesting interval
		if(($#interv_int == -1)or($other_ind[1] > $interv_int[$#interv_int][1]))
		{
		    $bool_local_not_recorded_interval = 1;
		}
		
		# for other seq, we have to find overrepresented seq to know wich position extend 
		($not_verbose)or print "On lance SORT_IND bool_first_call @bool_first_call sur les sequences indices @other_ind dans sub SORT avec border_words @$r_border_words et bool_local_not_recorded_interval $bool_local_not_recorded_interval, other_ind @other_ind, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice: line ".__LINE__."\n";

		&SORT_IND($other_ind[0], $other_ind[1], @bool_first_call , $bool_local_not_recorded_interval, $r_border_words, $r_proba_ext, 0, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice);
	    }
	    
	}
    }
    # END LOOP TO TREAT EVERY LOW PROBA, if beg_loop == end_loop, we treat every proba at pos proba_min, NOT EVERY **************************

    if(! $not_verbose)
    {
	print ' ' x length($string_let_pos_p);
	print "FIN SUB SORT_IND @_\n";
    }

}

# ************************************************************************************************

# ************************************************************************************************
# Compute a score corresponding to the rapport of number of corresponding motifs found in intergenic area on number of motifs found on the whole genome

sub RAPPORT_CRITERIA(\$\$\$\$\$)
{
    my ($r_regexp_motif, $r_rapport_score, $r_nb_in_promot, $r_nb_in_whole_2_sens) = @_;

    # ($not_verbose)or print "SUB RAPPORT_CRITERIA @_\n";
    $not_verbose or print "sub RAPPORT_CRITERIA at ",scalar(localtime(time())),"\n";

    # print "parameters $$r_regexp_motif, $way_to_whole_genome\n";
    

    
    # ($not_verbose)or print "motif perl $$r_regexp_motif motif grappe $regexp_motif_grappe, motif grappe inverse $reverse_regexp_motif_grappe\n";
    # ******* SEARCH FOR MOTIF IN WHOLE SEQUENCE **********************************************

    
    # file created to store interesting seq asociated with a motif
    # print "sp $sp\n";
    # print "name_fic_sp $name_fic_sp[$sp]\n";
    # print "r_regexp_motif $$r_regexp_motif\n";
    my $name_mot   = 'motif'.$name_fic_sp[$sp].'_'.$$r_regexp_motif.".txt";
    my $name_fsor  = $MOTIF_DIR_bact_min_nb_of_waited_sigma_motifs.$name_mot;
    my $name_fsor0 = $MOTIF_DIR_bact_min_nb_of_waited_sigma_motifs.'motif'.$name_fic_sp[$sp].'_'.$$r_regexp_motif."_0.txt";

    if(-e $name_fsor0)
    {
	($not_verbose)or print "$name_fsor0 already exists, not interesting\n";
	$$r_rapport_score = 0;
	# die "$name_fsor0 already exists, not interesting ln 4132\n";
    }
    elsif(-e $name_fsor)
    {
	($not_verbose)or print "$name_fsor already exists, interesting, we open it\n";
	# if rapport has already been computed for these species and motif, we get it in corresponding existing file
	open(FSOR,"< $name_fsor")or die "$prog_tag [Error] Impossible to open $name_fsor:$!, line ".__LINE__."\n";
	my $prem_ln  = (<FSOR>);
	my $sec_ln   = (<FSOR>);
	my $third_ln = (<FSOR>);
	if((! defined $prem_ln)or(! defined $prem_ln))
	{
	    die "$prog_tag [Error] prem_ln or $sec_ln not defined for $name_fsor file:$!, line ".__LINE__."\n";
	}
	close(FSOR)or die "$prog_tag [Error] Impossible to close FSOR:$!, line ".__LINE__."\n";
	$prem_ln =~ /rapport ([\d\.]+)/;
	$$r_rapport_score = $1;
	($not_verbose)or print "RAPPORT READ: dol1 $1\n";
	$sec_ln =~ /numerator (\d+)/;
	$$r_nb_in_promot = $1;
	$third_ln  =~ /divisor (\d+)/;
	$$r_nb_in_whole_2_sens = $1;
	# die "$name_fsor already exists, interesting, we open it ln 4154\n";
    }    
    # if rapport has never been computed for this species and motif, we compute it
    else
    {
	# open(PWD,'pwd |')or die "Impossible to open pwd command: $!\n";
	($not_verbose)or print "no existing file, we compute to create one\nWe do ls ".$way_to_whole_genome.$deb_file_whole_genome.$name_fic_sp[$sp]."*\n";
	my $pid;

	# die "no existing file, we compute to create one\nWe do ls ".$way_to_whole_genome.'not_fasta_whole_seq0_total_'.$name_fic_sp[$sp]."* ln 4163\nName will be $name_fsor or $name_fsor0\n";

        # ZONE DE SEGMENTATION FAULT *****************************************
	# print `free`."\n";

	my $sequence = '';
	my $bool_recup_seq_match = 1;
	# $$r_regexp_motif =~ /^(\w+)(\\w\{\d+,\d+\})(\w+)$/;

	# CREATION OF REVERSE COMPLEMENT OF REGULAR EXPRESSION******
	my (@normal_letters,@pos_normal_letter);
	
	while($$r_regexp_motif =~ /([ACGTacgt]+)(?:[^\]CGATacgt]|$)/g)
	{
	    push @normal_letters,  $1;
	    push @pos_normal_letter, $-[0];
	}
	my (@w_d,@pos_w_d);

# while($$r_regexp_motif =~ /(\\w(?:\{[\d\,]+\})|\[\w+\])/g)
	while($$r_regexp_motif =~ /(\\w(\{[\d\,]+\})?|\[\w+\])/g)
	{
	    my $el = $1;
	    $el =~ tr/ACGTacgt/TGCAtgca/;
	    push @w_d, $el;
	    push @pos_w_d, $-[0];
	}
	
	my $reverse_regexp_motif = '';
	
	my $ind_normal_letters = $#pos_normal_letter;
	my $ind_w_d = $#pos_w_d;
	my $bool_not_end_reverse = 1;
	do
	{
	    ($not_verbose)or print "pos_normal_letter de $ind_normal_letters: $pos_normal_letter[$ind_normal_letters]\npos_w_d de $ind_w_d: $pos_w_d[$ind_w_d]\n";
	    if(($ind_normal_letters>-1) and(($ind_w_d == -1)or($pos_normal_letter[$ind_normal_letters] > $pos_w_d[$ind_w_d])))
	    {
		($not_verbose)or print "on met ".&complement_seq($normal_letters[$ind_normal_letters])." ind $ind_normal_letters dans reverse_regexp_motif\n";
		$reverse_regexp_motif .= &complement_seq($normal_letters[$ind_normal_letters--]);
	    }
	    elsif(($ind_w_d>-1)and(($ind_normal_letters == -1)or($pos_normal_letter[$ind_normal_letters] < $pos_w_d[$ind_w_d])))
	    {
		($not_verbose)or print "on met $w_d[$ind_w_d] ind $ind_w_d dans reverse_regexp_motif\n";
		$reverse_regexp_motif .= $w_d[$ind_w_d--];
	    }
	    else{
		print "tab pos_normal_letter de taille $#pos_normal_letter: @pos_normal_letter\n********\n";
		print "ind_normal_letters $ind_normal_letters, ind_w_d $ind_w_d\n";
		print "pos_w_d de $ind_w_d vaut $pos_w_d[$ind_w_d]\n";
		print "pos_normal_letter de $ind_normal_letters vaut $pos_normal_letter[$ind_normal_letters]\n";
		system("echo \"\" | mail -s \"PB: $$r_regexp_motif bug in creating reverse complement\" $author_mail");
		exit;
	    }
	}until(($ind_normal_letters == -1)and($ind_w_d == -1));
	
	print "traitement motif $$r_regexp_motif, reverse_regexp $reverse_regexp_motif\n";
	# END CREATION OF REVERSE COMPLEMENT OF REGULAR EXPRESSION******


	# if(defined $3)
	# {
	    # ($not_verbose)or print "regexp $$r_regexp_motif dol1 $1 dol2 $2 dol3 $3\n";
	    # $reverse_regexp_motif = &complement_seq($3).$2.&complement_seq($1);
	# }
	# else
	# {
	  #  $$r_regexp_motif =~ /^(\w+)$/;
	  #  print "dol3 not defined, $$r_regexp_motif is one word\n";
	  #  $reverse_regexp_motif = &complement_seq($1);
	  #  if(!defined $1){ die "dol1 not defined for regexp_motif $$r_regexp_motif:$!"; }
	# }
	($not_verbose)or print "r_regexp_motif $$r_regexp_motif reverse_regexp $reverse_regexp_motif\n";
	 # die;

	# ****************************************************************************
	# BEGIN COUNT in WHOLESEQ ****************************************************
	# ****************************************************************************

	my $divisor = 0;

	# print "RAPPORT_CRITERIA\n";

	# my @files_grappe = (<LSWHOLE>);

	# print "LSWHOLE @files_grappe"."_toto\n";
	# $DB::single;
	 # while( my $f = <LSWHOLE>)
	{
	    my $res_l = 0;

	    # count whole seq sens 0
	    # open(F,$f)or die "Impossible to open $f file:$!";

	    # we go to beginning of whole genome file
	    seek(FWG,0,0);

	    # one line = complete genome OR plasmid, NOT a segmented sequence
	    while( my $whole_seq = <FWG>)
	    {
		# g = every occurrences, i = upper or lowercase
		while($whole_seq =~ /$$r_regexp_motif/gi)
		{
		    $res_l++;
		    if($res_l == $max_nb_of_waited_sigma_motifs)
		    {
			($not_verbose)or print "dans la boucle ligne ".__LINE__.", on trouve $$r_regexp_motif $res_l fois, TROP\n"; 
			$divisor = 0;
			# close F or die "Impossible to close F:$!";
			goto END_WHOLE_COUNT;
			# last COUNT1;
		    }
		}
	    }

	    my $added_count = 0;
	    
	    $divisor += $res_l; 
	    $added_count = $res_l;

	    # FIN DE ZONE DE SEGMENTATION FAULT *****************************************
	    ($not_verbose)or print "on ajoute res_l $res_l donne divisor1 $divisor\n"; 

	    if($$r_regexp_motif ne $reverse_regexp_motif)
	    {
		# count whole seq sens 1

		$res_l = 0;
		# we go to beginning of whole genome file
		seek(FWG,0,0);
		while( my $whole_seq = <FWG>)
		{
		    while($whole_seq =~ /$reverse_regexp_motif/gi)
		    {
			$res_l++;
			if($res_l == $max_nb_of_waited_sigma_motifs)
			{
			    ($not_verbose)or print "dans la boucle ligne ".__LINE__.", on trouve $reverse_regexp_motif $res_l fois, TROP\n"; 
			    $divisor = 0;
			   # close F or die "Impossible to close F:$!";
			    goto END_WHOLE_COUNT;
			   # last COUNT;
			}
		    }
		   
		}
		
		$divisor += $res_l; 
		($not_verbose)or print "divisor2 $divisor res_l $res_l\n";

	    }
	    else
	    {
		$divisor += $added_count; 
		($not_verbose)or print "motif PALYNDROMIQUE, on ajoute $added_count pour l'autre sens, donne $divisor\n";
	    }
	    if($res_l < $min_nb_of_waited_sigma_motifs)
	    {
		($not_verbose)or print "dans la boucle ligne ".__LINE__.", on trouve $$r_regexp_motif $res_l fois, TROP PEU\n"; 
		$divisor = 0;
	    }
	   # close F or die "Impossible to close F:$!";
	  END_WHOLE_COUNT:
	} # end COUNT
	# close LSWHOLE or die "Impossible to close LSWHOLE:$!";

	# ****************************************************************************
	# END COUNT in WHOLESEQ   ****************************************************
	# ****************************************************************************


	# ******* SEARCH FOR MOTIF IN WHOLE SEQUENCE **********************************************

	($not_verbose)or print "divisor $divisor\n";
	#    ($divisor)or print "on retourne zero car diviseur vaut 0\n";
	if(($divisor == 0)or($divisor > $max_nb_of_waited_sigma_motifs))
	{
	    $$r_rapport_score = 0;
	    # we have to record bad results not to compute it again
	    open(FSOR,">$name_fsor0")or die "$prog_tag [Error] Impossible to open $name_fsor0: $!, line ".__LINE__."\n";
	    print FSOR "rapport $$r_rapport_score\n";
	    print FSOR "numerator 0\n";
	    print FSOR "divisor 1\n";
	    close FSOR or die "$prog_tag [Error] Impossible to close FSOR:$!, line ".__LINE__."\n";
	    ($not_verbose)or print "rapports recorded: 0 in rapport_criteria\n";
	}
	else
	{
	    # we go to begin of seq files where we search for motif occurrences
	    seek($handle_ref_glob_sp0, 0, 0);
	    seek($handle_ref_glob_sp1, 0, 0);
    
	    ($not_verbose)or print "we go to beginning of $DATA_DIR".'SEQ_'.$name_fic_sp[$sp].".txt file pour comptage numerateur\n";
	    
	    # count upstr_seq
	    
	    my $numerator  = 0;
	    my $pos        = '';
	    my $real_motif = '';
	    my $local_cpt  = 0;
	    my @fsor       = (); #  tmp table to store what has to be written in $name_fsor file 

	    my @geneline = (); # stores every lines owning motif

	    my $nr_line = 0;

	    while(my $fseq = <$handle_ref_glob_sp0>)
	    {	
		# my $nr_line = $.;
		$nr_line++;
		my $fid = <$handle_ref_glob_sp1>;
		# print "ligne $nr_line four $fid, $fseq";

		# AJOUTER SI ($fid eq '') next; les sequences inexistantes correspondent aussi a des annotations, des lignes vides permettent d'avoir une equivalence
		# du nombre de lignes entre fichier de seq et annotations

		# my $tampon2 = '';
		# seek(FANNOT.$sp, $annot_byte_position[$sp][$nr_line], 0);
		# read FANNOT.$sp, $tampon2, $annot_byte_position[$sp][$nr_line+1] - $annot_byte_position[$sp][$nr_line];
		# print "$fid ligne $nr_line annot\n$tampon2\npos depart $annot_byte_position[$sp][$nr_line]\n";
		# exit;
		chomp($fseq);
		my $length_seq = length($fseq);
		my $bool_found = 0;
		
		while( $fseq =~ /($$r_regexp_motif)/gi )
		{  
		    $real_motif .= $1.' ';
		    # $pos .= ' -'.($length_seq - pos($fseq) + length($1)); # chgmt
		    $pos .= ' -'.($length_seq - $-[0]);

		    $numerator++;
		    $local_cpt++;
		    $bool_found = 1;

		    # print "on a le(s) motif(s) $real_motif numerator vaut $numerator, compteur local vaut $local_cpt\n";
		}
		($bool_found) and push @geneline, $nr_line;

		if(($bool_recup_seq_match)and($pos ne ''))
		{
		    # $fid =~ /:\s*(.+)$/;
		    # update 2020 07 21 to add chr/contig info, modif 2020 08 19 adding ^>
		    $fid =~ /^>(\S+):\s*(.+)$/;
		    # print "fid $fid";
		    # print "dol1 $1\n\n";
		    # push @fsor, ">$1 motif $real_motif nb motif $local_cpt, pos deb trad $pos\n$fseq\n";
		    # update 2020 07 21 to add chr/contig info
		    push @fsor, ">$1:$2 motif $real_motif nb motif $local_cpt, pos deb trad $pos\n$fseq\n";
		    # exit;
		    ($not_verbose)or print "pos vaut $pos, numerator $numerator pour $1\n";
		}
		if($pos ne '')
		{
		    # print "pos $pos pour $fseq\n";
		    $pos = '';
		}
		$real_motif = '';
		$local_cpt = 0;
		
	    }
	    # exit;
	    ($not_verbose)or print "numerator $numerator\n";
	    
	    # close FSORSEQ.$sp or die "Impossible to close FSORSEQ.$sp:$!";
	    # close FSORID.$sp or die "Impossible to close FSORID.$sp:$!";
	    # print "geneline @geneline\n";
	    # die "geneline @geneline\n";
         if($numerator > $divisor)	
	    {
print "PB: numerator $numerator greater than divisor $divisor\n";
		$divisor=$numerator;	
            }
else{
print "OK: numerator $numerator, $divisor\n";
}
        
	    $$r_rapport_score = $numerator/$divisor;
	 
	    # we edit interesting seq  and motifs only if they are found interesting only with orthologs according to criteria
	    if( $$r_rapport_score >= $rapports[$sp] )
	    {
		open(FSOR,">$name_fsor")or die "$prog_tag [Error] Impossible to open $name_fsor line ".__LINE__.": $!, line ".__LINE__."\n";
		print FSOR "rapport $$r_rapport_score\n";
		print FSOR "numerator $numerator\n";
		print FSOR "divisor $divisor\n";
		($bool_recup_seq_match)and print FSOR @fsor;
		close FSOR or die "$prog_tag [Error] Impossible to close FSOR:$!, line ".__LINE__."\n";

		# we open file where will be stored annotations of genes having motif in upstream region
		open(FANNOTMOT,">".$way_to_annot_motif_dir.'annot_'.$name_mot)or die "$prog_tag [Error] Impossible to open ".$way_to_annot_motif_dir."annot_$name_mot line ".__LINE__.": $!\n";
		my $tampon = '';

		# NOTE: EACH seq is on a line (without previous comment), so each box of geneline table
		# corresponds to the sequence
		# for last annot, we have recorded last byte (so, we have one byte value to avoid specific treatment for last case)

		for my $l (@geneline){
		   # print "boucle ".__LINE__."\n";

		    seek(FANNOT.$sp, $annot_byte_position[$sp][ $l-1 ], 0);
		    if($annot_byte_position[$sp][ $l-1 ] > $annot_byte_position[$sp][ $l ]){
			&verif_annot_byte_position($annot_byte_position[$sp]);
			die "$prog_tag [Error] for $name_mot annot, negative length, l vaut $l:  $annot_byte_position[$sp][ $l ] -  $annot_byte_position[$sp][ $l-1 ] vaut ".( $annot_byte_position[$sp][ $l ]- $annot_byte_position[$sp][ $l-1 ]).", line ".__LINE__."\n";
		    }
		    read FANNOT.$sp, $tampon, $annot_byte_position[$sp][ $l ] - $annot_byte_position[$sp][ $l-1 ];
		   # print "on lit la ligne $l, $tampon qui va de $annot_byte_position[$sp][ $l ] et fait en octet ".($annot_byte_position[$sp][ $l ] - $annot_byte_position[$sp][ $l-1 ])."\n********\n";
		    print FANNOTMOT $tampon,"\n";
		}

		# ouvrir et fermer fichier une seule fois
		# $index_annot_f[$i]
		close FANNOTMOT;

		$not_verbose or print "Le fichier genere s'appelle $name_fsor et comporte $numerator motifs trouves rapports $$r_rapport_score\n";
	    }
	    else
	    {
		open(FSOR,">$name_fsor0")or die "$prog_tag [Error] Impossible to open $name_fsor0: $!, line ".__LINE__."\n";
		print FSOR "rapport $$r_rapport_score\n";
		print FSOR "numerator $numerator\n";
		print FSOR "divisor $divisor\n";
		close FSOR or die "$prog_tag [Error] Impossible to close FSOR:$!, line ".__LINE__."\n";
		($not_verbose)or print "rapports $$r_rapport_score ds rapport_criteria\n";
	    }
	    # print "motif $$r_regexp_motif\nmotif inverse $reverse_regexp_motif\nmotif grappe $regexp_motif_grappe\nmotif inverse grappe $reverse_regexp_motif_grappe\nnumerator $numerator divisor $divisor\n";
	    
### THIS PART FOR THE MOMENT IS BEING PATCHED UP FOR DISCREPENCY IN VALUE of numerator and denom. LATER TFABRICE HAS TO UPDATE WITH NEW ALGORITHM
# 	    if($numerator > $divisor)
# 	    {
# 		$divisor=$numerator;		
# 		print "motif posant probleme dans PBmotif.txt\n";
# 		open(F,"> PBmotif.txt")or die "Impossible to open PBmotif.txt: $!";
# 		print F @fsor;
# 		close F or die "Impossible to close F:$!";
# 		die "PB numerator $numerator sup a divisor $divisor regexp $$r_regexp_motif: $!";
#	    }
	    $$r_nb_in_promot = $numerator;
	    $$r_nb_in_whole_2_sens = $divisor;

	}
    }
    # if($$r_regexp_motif =~ /ggaat/){ exit; }
}    

# ************************************************************************************************    
# ************************************************************************************************    

# create a regexp according to low proba < 0.2 found at 1 position
sub MATRICE_CREATE2(\@\@\$\@\@$$\$$){

    $not_verbose or print "sub MATRICE_CREATE2 @_\n";
    my $threshold_percent = 0.4;
   
    my ($r_border_words, $r_low_proba_pos_proba_min, $r_regexp_motif, $r_tempo_spacer_bounds, $r_global_cpt_score_pos_proba_min, $nb_seq, $nb_of_intern_l_added_to_boxes, $r_nb_l_added_matrice, $pos_proba_min) = @_;
    # print "r_nb_l_added_matrice BEG MATRICE $$r_nb_l_added_matrice, pos_proba_min $pos_proba_min\n";
    # variable permettant d'viter que les extensions puissent se chevaucher, voire empieter sur le trinuleotide adjacent
    my $extension_interne = 0;
    
    # we add seeds used for grouping of seq
    my $motif35 = $r_border_words->[0];
    my $motif10 = $r_border_words->[1];
    my $crochet = '';
    # print "border_words $r_border_words->[0], $r_border_words->[1], r_tempo_spacer_bounds $r_tempo_spacer_bounds->[0], $r_tempo_spacer_bounds->[1]\n";
    if($r_tempo_spacer_bounds->[1] - $r_tempo_spacer_bounds->[0] < 0){
	die "$prog_tag [Error] PB line ".__LINE__." with spacer_bounds, lower $r_tempo_spacer_bounds->[0], higher $r_tempo_spacer_bounds->[1]\n";
    }
    elsif($r_tempo_spacer_bounds->[1] - $r_tempo_spacer_bounds->[0]  > 10){
	print "CAUTION: in sub MATRICE_CREATE2, line ".__LINE__." with spacer_bounds, lower $r_tempo_spacer_bounds->[0], higher $r_tempo_spacer_bounds->[1]\n";
    }
    
    # foreach position
    # for (my $i = 0; $i <= 3; $i++)
    # {
    # if intern extension == lower spacer bound, we can not extend motif with a letter between 
    if(($extension_interne == $r_tempo_spacer_bounds->[0])and($pos_proba_min == 1 or $pos_proba_min == 2)){ die "$prog_tag [Error] How to treat boxes without spacer line ".__LINE__.": indices must be read according to low proba\n"; }
    
    # (defined $r_low_proba->[$i])or next;
    
    # push @{$low_proba[$pos]}, [$ind_l, $proba];
    # we have only one letter with a low proba at this position
    if($#$r_low_proba_pos_proba_min == 0){
	
	# print "indice de la lettre $r_low_proba_pos_proba_min->[0][0]\n";
	# print "nb seq for this lettre $#{ $r_global_cpt_score_pos_proba_min->[ $r_low_proba_pos_proba_min->[0][0] ] }\n";
	# print "threshold_percent $threshold_percent\n";
	if( $#{ $r_global_cpt_score_pos_proba_min->[ $r_low_proba_pos_proba_min->[0][0] ] }  <= $threshold_percent * $nb_seq){ return 0; }
	# else{ print "nb seq involved in proba $#{ $r_global_cpt_score_pos_proba_min->[ $r_low_proba_pos_proba_min->[0][0] ] } > $threshold_percent * $nb_seq, we add letter $rev_H_l[ $r_low_proba_pos_proba_min->[0][0] ]\n"; }
	# exit;
	if   ($pos_proba_min == 0){ $motif35 = $rev_H_l[ $r_low_proba_pos_proba_min->[0][0] ].$motif35; }
	elsif($pos_proba_min == 1){ $motif35 .= $rev_H_l[ $r_low_proba_pos_proba_min->[0][0] ]; $extension_interne++; }
	elsif($pos_proba_min == 2){ $motif10 = $rev_H_l[ $r_low_proba_pos_proba_min->[0][0] ].$motif10; $extension_interne++; }
	elsif($pos_proba_min == 3){ $motif10 .= $rev_H_l[ $r_low_proba_pos_proba_min->[0][0] ]; }
	
	$$r_nb_l_added_matrice++;
	# print "une lettre, 35 $motif35, 10 $motif10\n";
    }
    else{
	
	my $bool_found_l = 0;	    
	
	# foreach letter at this pos
	for (my $j = 0; $j <= $#$r_low_proba_pos_proba_min; $j++){
	    
	    if( $#{ $r_global_cpt_score_pos_proba_min->[ $r_low_proba_pos_proba_min->[$j][0] ] } <= $threshold_percent * $nb_seq){ next; }
	    # else{ print "nb seq involved in proba $#{ $r_global_cpt_score_pos_proba_min->[ $r_low_proba_pos_proba_min->[$j][0] ] } > $threshold_percent * $nb_seq, we add letter $rev_H_l[ $r_low_proba_pos_proba_min->[$j][0] ]\n"; }
	    
	    if(not $bool_found_l){
		if   ($pos_proba_min == 0){ $motif35 = ']'.$motif35; }
		elsif($pos_proba_min == 1){ $motif35 .= '['; }
		elsif($pos_proba_min == 2){ $motif10 = ']'.$motif10; }
		elsif($pos_proba_min == 3){ $motif10 .= '['; }
		$bool_found_l = 1;
	    }
	    # liste de couples lettre, proba pour cette position
	    # print "i $pos_proba_min, j $j,";
	    # print " ref $r_low_proba_pos_proba_min->[$j], ";
	    # print " val @{$r_low_proba_pos_proba_min->[$j]}\n";
	    if   ($pos_proba_min == 0){ $motif35 = $rev_H_l[ $r_low_proba_pos_proba_min->[$j][0] ].$motif35; }
	    elsif($pos_proba_min == 1){ $motif35 .= $rev_H_l[ $r_low_proba_pos_proba_min->[$j][0] ]; }
	    elsif($pos_proba_min == 2){ $motif10 = $rev_H_l[ $r_low_proba_pos_proba_min->[$j][0] ].$motif10; }
	    elsif($pos_proba_min == 3){ $motif10 .= $rev_H_l[ $r_low_proba_pos_proba_min->[$j][0] ]; }
	    
	}
	if($bool_found_l){
	    if($pos_proba_min == 0){ $motif35 = '['.$motif35; }
	    elsif($pos_proba_min == 1){ $motif35 .= ']'; $motif35 =~ s/\[(\w)\]/$1/; $extension_interne++; }
	    elsif($pos_proba_min == 2){ $motif10 = '['.$motif10; $extension_interne++; }
	    elsif($pos_proba_min == 3){ $motif10 .= ']'; $motif10 =~ s/\[(\w)\]/$1/; }
	    $$r_nb_l_added_matrice++;
	}
	
	     # print "plusieurs lettres, 35 $motif35, 10 $motif10\n";
    }    
	     
    my $space = '';

    # we test if 
    if($extension_interne == $r_tempo_spacer_bounds->[0]){ 
	($extension_interne != $r_tempo_spacer_bounds->[0])and $space = '\w{'.($r_tempo_spacer_bounds->[0] - $extension_interne - $nb_of_intern_l_added_to_boxes).'}'; }
    else{ $space = '\w{'.($r_tempo_spacer_bounds->[0] - $extension_interne - $nb_of_intern_l_added_to_boxes ).','.($r_tempo_spacer_bounds->[1] - $extension_interne - $nb_of_intern_l_added_to_boxes).'}'; }
    
    $$r_regexp_motif = $motif35.$space.$motif10;
    
    &REMPLACE_N_PAR_ANTISLASHW($r_regexp_motif);
    # avoid to have at[a]taaa motif for instance (gives atataaa)
    $$r_regexp_motif =~ s/\[(\w)\]/$1/g; 

    $not_verbose or print "SUB MATRICE2: r_regexp_motif $$r_regexp_motif\nboites $r_border_words->[0], $r_border_words->[1]\nesapces initiaux $r_tempo_spacer_bounds->[0], $r_tempo_spacer_bounds->[1], extension_interne $extension_interne, spacer initiaux: $r_tempo_spacer_bounds->[0] $r_tempo_spacer_bounds->[1], r_nb_l_added_matrice $$r_nb_l_added_matrice\n";

}



# ************************************************************************************************   
# ************************************************************************************************    

# create a regexp according to low proba < 0.2 found at 1 position
sub MATRICE_CREATE2_1L_AT_1POS(\@\@\$\@$\$$){

    $not_verbose or print "sub MATRICE_CREATE2_1L_AT_1POS @_\n";
   
    my ($r_border_words, $r_low_proba_pos_proba_min_1l, $r_regexp_motif, $r_tempo_spacer_bounds, $nb_of_intern_l_added_to_boxes, $r_nb_l_added_matrice, $pos_proba_min) = @_;
    # print "r_nb_l_added_matrice BEG MATRICE $$r_nb_l_added_matrice, pos_proba_min $pos_proba_min\n";
    # variable permettant d'viter que les extensions puissent se chevaucher, voire empieter sur le trinuleotide adjacent
    my $extension_interne = 0;
    
    # we add seeds used for grouping of seq
    my $motif35 = $r_border_words->[0];
    my $motif10 = $r_border_words->[1];

    # print "border_words $r_border_words->[0], $r_border_words->[1], r_tempo_spacer_bounds $r_tempo_spacer_bounds->[0], $r_tempo_spacer_bounds->[1]\n";
    if($r_tempo_spacer_bounds->[1] - $r_tempo_spacer_bounds->[0] < 0){
	die "$prog_tag [Error] PB line ".__LINE__." with spacer_bounds, lower $r_tempo_spacer_bounds->[0], higher $r_tempo_spacer_bounds->[1]\n";
    }
    elsif($r_tempo_spacer_bounds->[1] - $r_tempo_spacer_bounds->[0]  > 10){
	print "CAUTION: in sub MATRICE_CREATE2_1L_AT_1POS, line ".__LINE__." with spacer_bounds, lower $r_tempo_spacer_bounds->[0], higher $r_tempo_spacer_bounds->[1]\n";
    }
    
    # if intern extension == lower spacer bound, we can not extend motif with a letter between 
    if(($extension_interne == $r_tempo_spacer_bounds->[0])and($pos_proba_min == 1 or $pos_proba_min == 2)){ die "$prog_tag [Error] How to treat boxes without spacer line ".__LINE__.": indices must be read according to low proba\n"; }
    
    # print "indice de la lettre $r_low_proba_pos_proba_min_1l->[0]\n";
    # print "nb seq for this lettre $#{ $r_global_cpt_score_pos_proba_min->[ $r_low_proba_pos_proba_min_1l->[0] ] }\n";
    # print "threshold_percent $threshold_percent\n";

    if   ($pos_proba_min == 0){ $motif35 = $rev_H_l[ $r_low_proba_pos_proba_min_1l->[0] ].$motif35; }
    elsif($pos_proba_min == 1){ $motif35 .= $rev_H_l[ $r_low_proba_pos_proba_min_1l->[0] ]; $extension_interne++; }
    elsif($pos_proba_min == 2){ $motif10 = $rev_H_l[ $r_low_proba_pos_proba_min_1l->[0] ].$motif10; $extension_interne++; }
    elsif($pos_proba_min == 3){ $motif10 .= $rev_H_l[ $r_low_proba_pos_proba_min_1l->[0] ]; }
	
    $$r_nb_l_added_matrice++;

    # print "une lettre, 35 $motif35, 10 $motif10\n";
   
    my $space = '';

    # we test if 
    if($extension_interne == $r_tempo_spacer_bounds->[0]){ 
	($extension_interne != $r_tempo_spacer_bounds->[0])and $space = '\w{'.($r_tempo_spacer_bounds->[0] - $extension_interne - $nb_of_intern_l_added_to_boxes).'}'; }
    else{ $space = '\w{'.($r_tempo_spacer_bounds->[0] - $extension_interne - $nb_of_intern_l_added_to_boxes ).','.($r_tempo_spacer_bounds->[1] - $extension_interne - $nb_of_intern_l_added_to_boxes).'}'; }
    
    $$r_regexp_motif = $motif35.$space.$motif10;
    
    &REMPLACE_N_PAR_ANTISLASHW($r_regexp_motif);

    # print "SUB MATRICE_CREATE2_1L_AT_1POS: r_regexp_motif $$r_regexp_motif\nboites $r_border_words->[0], $r_border_words->[1]\nesapces initiaux $r_tempo_spacer_bounds->[0], $r_tempo_spacer_bounds->[1], extension_interne $extension_interne, spacer initiaux: $r_tempo_spacer_bounds->[0] $r_tempo_spacer_bounds->[1], r_nb_l_added_matrice $$r_nb_l_added_matrice\n";
}



# ************************************************************************************************ 

 
# ************************************************************************************************

# create a matrice based on seq related to current interval, create a position specific matrice and deduct regular expression motif wich will be used for ratio computing
sub MATRICE_CREATE($$\@\$\$\@)
{
    my ($last_ind_seq_set_first, $last_ind_seq_set_sec, $r_tmp_display, $r_regexp_motif, $r_nb_l_added_matrice, $r_tempo_spacer_bounds) = @_;
    
    $not_verbose or print "sub MATRICE_CREATE at ",scalar(localtime(time())),"\n";
    # **************************************************************************************
    # A CE STADE, LES INDICES DES SEQUENCES SI ELLES ETAIENT ORDONNÉES SONT STOCKES DANS SORTED_IND, LES SEQUENCES EN
    # ELLES-MEMES NE SONT PAS ORDONNÉES POUR EVITER UNE PERTE DE TEMPS
    
    # print "SUB MATRICE_CREATE @_\n";
    # ($initial_border_words[1] =~ /ggta/)and print "ggta: SUB MATRICE_CREATE @_, sp $sp\n";
    # pseudo-stat on letters surrounding motif boxes
    # this first constant fixes the width of the number of letters we are taking in account for
    # pseudo-stat computing each side of each trinucleotid
    my $WIDTH_BORDER_STAT = 4;
    

    # TREATMENT AND DISPLAY FOR EACH INTERESTING INTERVAL 
    
    my $deb_prev_intit = ' ';    

    # treatment of our interval**************************************
    
    # tab_count ininitialization
    my @tab_count = ();
    for my $j(0..4*$WIDTH_BORDER_STAT-1)
    {
	@{$tab_count[$j]} = (0,0,0,0);
    }
    
    my $subs = ''; # tempo var for counted substring
    my @tab_correctif_somme; # the seq is to short, no lettre corresponds, so total number of computed lettres will be corrected (in comparison with number of sequences)(first part corresponds to letters on the left of -35 box, second part corresponds to letters on the right of -10 box) 

    # CORRECTION OF SUMMS IN CASE OF TO SHORT SEQ***********************************************
    # if($sp){ $last_ID = $#{$seq_sp[$sp][1]}; }
    # else   { $last_ID = $#{$seq_sp[$sp][0]}; }
    
    @tab_correctif_somme = (0)x (2*$WIDTH_BORDER_STAT); # corrections
    
    # ($not_verbose)or print "VERIF matrice\n";
    for my $cpt_l($last_ind_seq_set_first..$last_ind_seq_set_sec)
    {
	my $sorted_i = $sorted_ind[$cpt_l];

	# ($not_verbose)or print "Seq concernee $seq_sp[$sp][ $sorted_i ] d'ind $cpt_l dans sorted_ind et d'ind reel $sorted_i\n";    
	for my $id_l_red(1..$WIDTH_BORDER_STAT)
	{
	    # print "cpt_l $cpt_l, sorted_i $sorted_i\n";
	    # print "seq_sp ".$seq_sp[$sp][ $sorted_i ]."\n";
	    # print "pos_tri $pos_tri[$sp][ $sorted_i ][0]\n";

	    # left borders of left trinuc
	    if( (defined ($subs = substr($seq_sp[$sp][ $sorted_i ],$pos_tri[$sp][ $sorted_i ][0]-$id_l_red,1)))and($subs ne '') )
	    {
		$tab_count[$WIDTH_BORDER_STAT-$id_l_red][$H_l{$subs}]++;
	    }
	    else
	    {
		$tab_correctif_somme[$WIDTH_BORDER_STAT-$id_l_red]++;
	    }
	    # left border of right trinuc
	    $tab_count[3*$WIDTH_BORDER_STAT-$id_l_red][$H_l{ substr($seq_sp[$sp][ $sorted_i ],$pos_tri[$sp][ $sorted_i ][1]-$id_l_red,1) }]++;
	    
	    # right borders of left trinuc
	    $tab_count[$WIDTH_BORDER_STAT+$id_l_red-1][$H_l{ substr($seq_sp[$sp][ $sorted_i ],$pos_tri[$sp][ $sorted_i ][0]+$LOCAL_WIDTH_OF_L_HBOXES+$id_l_red-1,1) }]++;

	    # position at right of right trinuc ??
	    my $ind_l_seq = $pos_tri[$sp][ $sorted_i ][1]+$LOCAL_WIDTH_OF_R_HBOXES+$id_l_red-1;

	    # right border of right trinuc
	    if( ( length($seq_sp[$sp][ $sorted_i ]) > $ind_l_seq )and(defined ($subs = substr($seq_sp[$sp][ $sorted_i ], $ind_l_seq, 1)))and($subs ne '') )
	    {
		$tab_count[3*$WIDTH_BORDER_STAT+$id_l_red-1][$H_l{$subs}]++;
	    }
	    else
	    {
		$tab_correctif_somme[$WIDTH_BORDER_STAT+$id_l_red-1]++;
	    }
	    
	}
    }
    # END CORRECTION OF SUMMS IN CASE OF TO SHORT SEQ***********************************************
    
    # COMPUTE OF RATES DISPLAY**********************************************************************
    
    my $nb_seq = $last_ind_seq_set_sec - $last_ind_seq_set_first + 1; 
    
    @$r_tmp_display = ('a  :','c  :','g  :','t  :'); # table used for stat display
    my @stat = ();
    # my @motif_letters = ();	
    
    die "$prog_tag [Error] PB line ".__LINE__."\n" if($nb_seq == 0);
    for my $j(0..4*$WIDTH_BORDER_STAT-1)
    {
	# print "i $i j $j acgt @{$tab_count[$i][$j]}\n";
	for my $k(0..$#$r_tmp_display) # lines of DISPLAY correspond to letters
	{
	    # print "seq $seq_sp[$sp][$i][$sorted_i]\n";
	    if($j == $WIDTH_BORDER_STAT)
	    {
		# ??CHANGER
		$r_tmp_display->[$k] .= $initial_border_words[0].' ';
	    }
	    elsif( $j == ($WIDTH_BORDER_STAT << 1))
	    {
		$r_tmp_display->[$k] .= "   ";
	    }
	    elsif( $j == $WIDTH_BORDER_STAT*3)
	    {
		# ??CHANGER
		$r_tmp_display->[$k] .= $initial_border_words[1].' ';
	    }
	    
	    my $divisor;
	    if($j < $WIDTH_BORDER_STAT)
	    { 
		$divisor = $nb_seq -$tab_correctif_somme[$j]; 
	    }
	    elsif($j >= 3*$WIDTH_BORDER_STAT)
	    { 
		$divisor = $nb_seq -$tab_correctif_somme[$j-2*$WIDTH_BORDER_STAT]; 
	    }
	    else
	    { 
		$divisor = $nb_seq; 
	    }
	    
	    if($divisor > 0)
	    {
		$stat[$j][$k] = $tab_count[$j][$k]/$divisor;
		# $motif_letters[$j] = ;
		$r_tmp_display->[$k] .= sprintf "%2.2f ",$stat[$j][$k];
	    }
	}
	
    }
    # END COMPUTE OF RATES DISPLAY**********************************************************************
    
    # NB_L_ADDED_MATRICE COMPUTING**********************************************************************************
    
    my $motif35_to_search = '';
    my $motif10_to_search = '';
    
    # $intitules_ind_prem_l =~ /subsp (\d+)(\s+.*?(\d+)?)?$/;
    
    # first case countains lower limit for spacer between the two boxes
    # second case countains higher limit for spacer between the two boxes
    # my @spacer_bounds = ();
    
    # if(! defined $1){ die "dol1 non defini pour $intitules_ind_prem_l\n"; }
    
    # if( defined $3){ push @spacer_bounds, ($1,$3); ($not_verbose)or print "dol3 defini\n";}
    # else{            push @spacer_bounds, ($1,$1); ($not_verbose)or print "dol3 NON defini\n";}
    
    ($not_verbose)or print "INITIAL spacer_bounds @$r_tempo_spacer_bounds\n";

    # print "intitule $intitules_ind_prem_l dol1 $1 dol2 $2\n";

    # variable permettant d'viter que les extensions puissent se chevaucher, voire empieter sur le trinuleotide adjacent
    my $limit_extension_interne = $r_tempo_spacer_bounds->[0];

    # before -35 box

    # parameter COUNT_SURROUNDING_RATES(array of stat, 0, $WIDTH_BORDER_STAT-1, boolean to know if extension is from left to right (0) or from right to left (1 (reverse mode)), threshold for extension $r_nb_l_added_matrice, \@tab_count, ref to maximal extension possible not to overlap neighboring box (and previsous enxtension), boolean to know if extension is between the two boxes (1) or not (0))

    # for the gapped trinuc, added letters won't be significatives (other else it is not usefull as we already found results with trinuc without gaps...), so, we don't take in account normally added letters. 
    if($initial_border_words[0] =~ /n/)
    {
	$motif35_to_search = $initial_border_words[0];
    }
    else
    {
	$motif35_to_search .= &COUNT_SURROUNDING_RATES(\@stat, 
						       0, 
						       $WIDTH_BORDER_STAT-1, 
						       1, 
						       $r_nb_l_added_matrice, 
						       \@tab_count, 
						       \$limit_extension_interne, 
						       0).$initial_border_words[0];
    }
    defined $motif35_to_search or die "$prog_tag [Error] motif35_to_search not defined line ".__LINE__."\n";
    
    # ($initial_border_words[1] =~ /ggta/)and print "r_nb_l_added_matrice $$r_nb_l_added_matrice after 35 left\n";
    # after 35 box
    if($limit_extension_interne > 0)
    {
	# only add to compense next substraction
	foreach(@$r_tempo_spacer_bounds){ $_ += length($motif35_to_search); }
	
	if($initial_border_words[0] !~ /n/)
	{
	    $motif35_to_search .= &COUNT_SURROUNDING_RATES(\@stat, 
							   $WIDTH_BORDER_STAT, 
							   2*$WIDTH_BORDER_STAT-1, 
							   0, 
							   $r_nb_l_added_matrice, 
							   \@tab_count, 
							   \$limit_extension_interne, 
							   1);
	}
	defined $motif35_to_search or die "$prog_tag [Error] motif35_to_search not defined line ".__LINE__."\n";
	
	# ($initial_border_words[1] =~ /ggta/)and print "r_nb_l_added_matrice $$r_nb_l_added_matrice after 35 right\n";
	foreach(@$r_tempo_spacer_bounds){ $_ -= length($motif35_to_search); }
	
	($not_verbose)or print "motif35_to_search $motif35_to_search subw1 $initial_border_words[0] spacer_bounds @$r_tempo_spacer_bounds\nFIN TRAITEMENT MOTIF 35 DEB TRAITEMENT MOTIF 10 nb_l_added_matrice $$r_nb_l_added_matrice\n";
	# print "motif35_to_search $motif35_to_search subw1 $initial_border_words[0] spacer_bounds @$r_tempo_spacer_bounds\nFIN TRAITEMENT MOTIF 35 DEB TRAITEMENT MOTIF 10 nb_l_added_matrice $$r_nb_l_added_matrice\n";
	
	# we add space to search motif between two seq
	
	# before 10 box
	($initial_border_words[1] =~ /n/)or $motif10_to_search .= &COUNT_SURROUNDING_RATES(\@stat, 2*$WIDTH_BORDER_STAT, 3*$WIDTH_BORDER_STAT-1, 1, $r_nb_l_added_matrice, \@tab_count, \$limit_extension_interne, 1);
	
	# ($initial_border_words[1] =~ /ggta/)and print "r_nb_l_added_matrice $$r_nb_l_added_matrice after 10 left\n";
	# ($not_verbose)or print "motif10_to_search $motif10_to_search juste apres nb_l_added_matrice $$r_nb_l_added_matrice\n";
	# print "motif10_to_search $motif10_to_search juste apres nb_l_added_matrice $$r_nb_l_added_matrice\n";
	
	if($motif10_to_search ne '')
	{
	    foreach(@$r_tempo_spacer_bounds){ $_ -= length($motif10_to_search); } 
	}
	
    }
    
    # after 10 box
    $motif10_to_search .= $initial_border_words[1];
    
    # ($not_verbose)or print "motif10_to_search $motif10_to_search apres AJOUT subw2\n";
    # print "motif10_to_search $motif10_to_search apres AJOUT subw2\n";
    
    ($initial_border_words[1] =~ /n/)or $motif10_to_search .= &COUNT_SURROUNDING_RATES(\@stat, 3*$WIDTH_BORDER_STAT, 4*$WIDTH_BORDER_STAT-1, 0, $r_nb_l_added_matrice, \@tab_count, \$limit_extension_interne, 0);
    
    # ($initial_border_words[1] =~ /ggta/)and print "r_nb_l_added_matrice $$r_nb_l_added_matrice after 10 right\n";
    # ($not_verbose)or print "subw2 $initial_border_words[1] motif10_to_search $motif10_to_search nb_l_added_matrice $$r_nb_l_added_matrice\n";
    
    # ($not_verbose)or print "SEARCHED WORD: limits @$r_tempo_spacer_bounds tirees de $intitules_ind_prem_l\n";
    if($r_tempo_spacer_bounds->[1] > 0)
    {
	$$r_regexp_motif = $motif35_to_search.'\w{'.$r_tempo_spacer_bounds->[0].','.$r_tempo_spacer_bounds->[1].'}'.$motif10_to_search;
	# ($initial_border_words[1] =~ /ggta/)and print "ggta: normal case for regexp generation\n";
    }
    else
    {
	# two boxes are only one in fact, without variable spacer
	$$r_regexp_motif = $motif35_to_search.$motif10_to_search;
	# ($initial_border_words[1] =~ /ggta/)and print "ggta: word fusion case for regexp generation\n";
    }
    
    &REMPLACE_N_PAR_ANTISLASHW($r_regexp_motif);
    # print "regexp_motif $$r_regexp_motif\n";
}




# complementSequence (seq) returns the complement of nucleotid sequence (seq). Nucleotids may be uppercase or lowercase.
# DOES NOT CHECK SEQUENCE FOR ILLEGAL NUCLEOTIDS.
# CAUTION, THERE IS A MAX SIZE for string passed as parameter
sub complementSequence($) 
{
    my $dnaSequence = shift @_;
    my $reverseSeq;
    
    $reverseSeq = reverse "$dnaSequence"; # Sequence is reversed
    $reverseSeq =~ tr/atgcATGC/tacgTACG/; # Individual nucleotids are complemented
    # ($not_verbose)or print "SUB Complement seq of $dnaSequence OK\n";
    return $reverseSeq;
}


# COMPUTE prob series **************************************************** 
# param(p,n) 
sub Cnk_p_q($$)
{
    my ($p,$n) = @_;
    my $n_prim = $n;
    my $q = 1-$p;
    my ($ref, $i, $j);
    my @puis_deux = (0)x ($n+1); 
    my @fct_gener = @puis_deux; 
    @puis_deux[0..1] = ($q,$p);
    $fct_gener[0] = 1.00;
    my @polynomeG = (); 
    my $r_pd = \@puis_deux;
    my $r_fg = \@fct_gener;
    my $r_pG = \@polynomeG;
    
    # ($not_verbose)or print "SUB Cnk_p_q @_\n";

    while($n_prim > 0)
    {
	if($n_prim & 1)
	{
	    # fct_gener = fct_gener*puis_deux
	    for $i(0..$n)
	    {   
		# print "i $i\n";
		$r_pG->[$i] = 0;

		for $j(0..$i){   $r_pG->[$i] += $r_pd->[$i-$j] * $r_fg->[$j]; }
	    }                                                                                                                                                     
	    $ref = $r_pG;
	    $r_pG = $r_fg;
	    $r_fg = $ref;
	}
	# die "@$r_fg\n";

	# puis_deux = puis_deux * puis_deux; 
	for $i(0..$n)
	{
	    # print "i $i\n";
	    $r_pG->[$i] = 0;
	    for $j(0..$i)
	    {
		$r_pG->[$i] += $r_pd->[$i-$j] * $r_pd->[$j];
	    }
	}
	$ref = $r_pG;
	$r_pG = $r_pd;
	$r_pd = $ref;
	$n_prim = $n_prim >> 1;
    }
    
    return $r_fg;
}
# ****************************************************************

# ****************************************************************
# param: 
# ref to proba table, 
# ref to trinuc (left extension of left trinuc) 
# ref to trinuc (right extension of left trinuc) 
# ref to trinuc (left extension of right trinuc) 
# ref to trinuc (right extension of right trinuc) 
# n for series computing
# 4 boolean, 1 by position (1 means extended, necessary update)
sub INIT_PROBA_TAB(\@\@$)
{
    my ($r_proba_ext, $r_border_words, $n) = @_;
    # my ($r_proba_ext, $r_border_words, $n, @pos_for_MAJ) = @_;
    
    $not_verbose or print "PROBA evalues on @$r_border_words, nbseq $n\n";
    # ($not_verbose)or print "SUB INIT_PROBA_TAB @_\n";
    
    $n++; # we had higher ind and we want cells number
    # ($not_verbose)or print "r_pos_for_MAJ @$r_pos_for_MAJ  a l'interieur de MATRICE\n";
    for my $pos(0..$#bool_first_call)
    {
	# ($pos_for_MAJ[ $pos ])or next;
	($bool_first_call[$pos]) or next;

	my $trinuc_proba; # trinuc used for proba compute (case when trinuc box is larger than three trinuc)
	# if pos mod 2 is 1, extension concerns letters on le right of trinuc box, we take # in ww###

	if   ($pos == 1){ $trinuc_proba = substr($r_border_words->[ 0 ], length($r_border_words->[ 0 ]) - 3, 3); }
	elsif($pos == 3){ $trinuc_proba = substr($r_border_words->[ 1 ], length($r_border_words->[ 1 ]) - 3, 3); }
	elsif($pos == 0){ $trinuc_proba = substr($r_border_words->[ 0 ], 0, 3); }
	elsif($pos == 2){ $trinuc_proba = substr($r_border_words->[ 1 ], 0, 3); }
	else            { die "$prog_tag [Error] Pos unknown for extension line ".__LINE__."\n"; }

	# if($pos % 2){ $trinuc_proba = substr($r_border_words->[ $pos ], 0, 3); }
	# if pos mod 2 is 0, extension concerns letters on le left of trinuc box, we take # in ###ww
	# else        { $trinuc_proba = substr($r_border_words->[ $pos ], length($r_border_words->[ $pos ]) - 3, 3); }

	for my $l(0..3)
	{
	    if(not defined $H_MM4_extension{ $trinuc_proba }[$sp][ $pos & 1 ][ $l ])
	    {
		die "$prog_tag [Error] H_MM4_extension of: trinuc_proba ($trinuc_proba), sp ($sp), pos ($pos) & 1, l ($l), not defined line ".__LINE__."\n";
	    }
	    if(not defined $n)
	    {
		die "$prog_tag [Error] number of sequences n not defined line ".__LINE__."\n";
	    }
	    $r_proba_ext->[ $pos ][ $l ] = &Cnk_p_q($H_MM4_extension{ $trinuc_proba }[$sp][ $pos & 1 ][ $l ], $n); # &Cnk_p_q(p,n) 
	    # ($not_verbose)or print "INTERIEUR @{ $r_proba_ext->[ $pos ][ $l ] }\n";
	}
    }

    if($#$r_proba_ext >= $n)
    {
	for my $pos(0..$#bool_first_call)
	{
	    for my $l(0..3)
	    {
		$#{ $r_proba_ext->[ $pos ][ $l ] } = $n-1;
		# ($not_verbose)or print "INTERIEUR2 @{ $r_proba_ext->[ $pos ][ $l ] }\n";
	    } 
	}
    }
} 

# ****************************************************************

# ****************************************************************
sub COUNT_NB_SEQ(\$\$$$$)
{

    my ($r_H_to_count_seq_by_spacer, $r_count_tmp, $pack_subsp_v, $subsp_v, $bool_ID_table) = @_;
   # ($not_verbose)or print "********************************************\nARGV COUNT_NB_SEQ @_\n";
    
    # print "table: we increase r_count_tmp whose value was $$r_count_tmp ";

    # count if table (EVENTUALLY number of boxes)
    if(($bool_ID_table)and(exists $r_H_to_count_seq_by_spacer->{$pack_subsp_v}))
    {
	$$r_count_tmp += scalar(@{$r_H_to_count_seq_by_spacer->{$pack_subsp_v}});
	($not_verbose)or print "and is then $$r_count_tmp (case 1)\n";
	# die "$!";
    }
    # count if H keys
    elsif(not $bool_ID_table)
    {
	while( my ($pack_num) = each %{$r_H_to_count_seq_by_spacer->{$pack_subsp_v}})
	{
	    ($not_verbose)or print "on traite ".unpack('N',$pack_num)." pour $subsp_v\n";
	    (exists $r_H_to_count_seq_by_spacer->{$pack_subsp_v}{$pack_num}[ 1 ])or next;
	    
	    if($r_H_to_count_seq_by_spacer->{$pack_subsp_v}{$pack_num}[ 1 ] =~ /\-/)
	    {
		die "$prog_tag [Error] - alors qu'il contient  $r_H_to_count_seq_by_spacer->{$pack_subsp_v}{$pack_num}[ 1 ]:$!, line ".__LINE__."\n"; # \nmotif $important_motifs[ unpack('N',$pack_num) ][0] ".unpack('U',$important_motifs[ unpack('N',$pack_num) ][1])." $important_motifs[ unpack('N',$pack_num) ][2]\n";
	    }
	    
	    for my $c( 0..length( $r_H_to_count_seq_by_spacer->{$pack_subsp_v}{$pack_num}[ 1 ] ) -1 )
	    {
		($not_verbose)or print "H_nbw: on ajoute ".substr( $r_H_to_count_seq_by_spacer->{$pack_subsp_v}{$pack_num}[ 1 ] , $c, 1)." extrait de ".$r_H_to_count_seq_by_spacer->{$pack_subsp_v}{$pack_num}[ 1 ]." a count_tmp $$r_count_tmp ";
		$$r_count_tmp += substr( $r_H_to_count_seq_by_spacer->{$pack_subsp_v}{$pack_num}[ 1 ] , $c, 1);
		($not_verbose)or print "qui vaut alors $$r_count_tmp\n";
	    }
	}
	($not_verbose)or print "and is then $$r_count_tmp (case 2)\n";
	# die "$!";
    }
    else
    {
	($not_verbose)or print "\n";
    }

    ($not_verbose)or print "END COUNT_NB_SEQ @_\n********************************************\n";
}
# ****************************************************************

sub TREATMENT_COUNT_NB_SEQ(\%$$$$)
{
    my ($r_H_to_count_seq_by_spacer, $count_tot_ind_tot, $num_subw1, $num_subw2, $pack_record_sp) = @_;

    # each motif ID 
    while( my ($pack_num) = each  %{$r_H_to_count_seq_by_spacer->{$pack_record_sp}})
    {
	# if at least one seq has been counted
	if( exists $r_H_to_count_seq_by_spacer->{$pack_record_sp}{$pack_num}[ 1 ] )
	{
	    # we record in H_nbw the data needed to browse all the records corresponding to trinuc1 trinuc2 spacer
	    $H_nbw{ $count_tot_ind_tot }{$num_subw1}{$num_subw2}{$pack_record_sp}{$pack_num} = $r_H_to_count_seq_by_spacer->{$pack_record_sp}{$pack_num}; 
	   # print "num ".unpack('U', $pack_num)."\npour pos sous-mots ".@{$r_H_to_count_seq_by_spacer->{$pack_record_sp}{$pack_num}[ 0 ]}.' et pour comptages '.$r_H_to_count_seq_by_spacer->{$pack_record_sp}{$pack_num}[ 1 ]."\n" ;
	}
	# else
	# {
	  #  die "r_H_to_count_seq_by_spacer does not exist\n";
	# }
    }
}

# ****************************************************************
sub TEST_NBSEQ_PER_INTERV
{
    # print "ARGV TEST_NBSEQ_PER_INTERV @_\n";
    my ($r_H_to_count_seq_by_spacer, $indice_middle_spacershift_table, $r_valid_int_sp, $r_H_nbw, $num_subw1, $num_subw2);

    my $bool_IDtable;
    # treatment if H_nbw: we count seq thanks scalar keys(%num)
    if($#_ == 4)
    {
	# (\%$\%$$)
	($r_H_to_count_seq_by_spacer, $indice_middle_spacershift_table, $r_H_nbw, $num_subw1, $num_subw2) = @_;
	$bool_IDtable = 0;
    }
    # treatment table: count is done according to number of recorded ids scalar(@id)
    else
    {
	# (\%$\%)
	($r_H_to_count_seq_by_spacer, $indice_middle_spacershift_table, $r_valid_int_sp) = @_;
	$bool_IDtable = 1;
    }

    # ($not_verbose)or print "subw2 $table_trinuc[$num_subw2] parcours r_trinuc1_trinuc2_spacer line ".__LINE__."..\n";
    my @count = ();
    my @count_tot = ();
    my $count_summ = 0;
    my $prev_subsp = '';
    my $read_keys_count = 0; # count met spacements (used only to know if we read the last key for spacement (to record in this case))
    my $tot_nb_keys = scalar(keys %{$r_H_to_count_seq_by_spacer}); # numbers of keys corresponding to spacement
    
    
    # read by increasing order of space width
  COUNTLOOP:for my $pack_subsp(sort unpack_numerique_U keys %{$r_H_to_count_seq_by_spacer})
  {
      my $subsp = unpack('U',$pack_subsp);
      
      ($not_verbose)or print "subsp $subsp parcours r_trinuc1_trinuc2_spacer line ".__LINE__."..\n";
      
      $read_keys_count++;
      
      if($#count == -1)
      {
	  $prev_subsp = $subsp;
	  $count_summ = 0;
	  
	  # initialization of @count
	  @count[0..$indice_middle_spacershift_table-1] = (0) x $indice_middle_spacershift_table;
	  
	  for my $subsp_v($subsp..$subsp + $indice_middle_spacershift_table)
	  {
	      my $pack_subsp_v = pack('U',$subsp_v);

	      # only for case r_H_to_count_seq_by_spacer is an H_table
	      if((not $bool_IDtable)and( scalar( keys( %{$r_H_to_count_seq_by_spacer->{$pack_subsp_v}} )) == 0 ))
	      {
		  push @count, 0;
	      }
	      else
	      {
		  my $count_tmp = 0;
		  # VERIF ??
		  &COUNT_NB_SEQ($r_H_to_count_seq_by_spacer, \$count_tmp, $pack_subsp_v, $subsp_v, $bool_IDtable);

		  push @count, $count_tmp;
		  $count_summ += $count_tmp;
		  ($not_verbose)or print "count_summ vaut $count_summ for $subsp_v\n";
	      }
	  }
	  @count_tot = ($count_summ) x @count;
	  ($not_verbose)or print "Initialisation count_tot @count_tot count @count\n";
      } # endif($#count == -1)
      else
      {
	  # shift number for table
	  my $decallage = 0;
	  ($prev_subsp eq '')or $decallage = $subsp - $prev_subsp;
	  
	  # we have to record the count for trinuc which are going to be erase at this turn
	  my $end_b = $decallage - 1;
	  if(($prev_subsp ne '')and($#count <= $decallage))
	  {
	      $end_b = $#count; 
	      ($not_verbose)or print "end_b mis a $end_b\n"; 
	  }
	  
	  ($not_verbose)or print "decallage vaut $decallage end_b vaut $end_b\n";
	  
	  for my $ind_tot(0..$end_b)
	  {
	      if( ($count[$ind_tot] !=0)and($count_tot[$ind_tot] >= $min_nb_of_waited_sigma_motifs) )
	      {
		  my $pack_record_sp = pack('U',$prev_subsp - ($indice_middle_spacershift_table - $ind_tot) );
		  if($bool_IDtable)
		  {
		      #pack record sp!!!!!!
		      $r_valid_int_sp->{ $pack_record_sp } = 1;
		     # $r_valid_int_sp->{ $prev_subsp - ($indice_middle_spacershift_table - $ind_tot) } = 1;
		      ($not_verbose)or print "r_valid_int_sp de ".($prev_subsp - ($indice_middle_spacershift_table - $ind_tot))." cree\n";
		  }
		  else
		  {
		      &TREATMENT_COUNT_NB_SEQ($r_H_to_count_seq_by_spacer, $count_tot[$ind_tot], $num_subw1, $num_subw2, $pack_record_sp);
		      # print $r_H_nbw->{ $count_tot_ind_tot }{$num_subw1}{$num_subw2}{$pack_record_sp}{$pack_num} = $r_H_to_count_seq_by_spacer->{$pack_record_sp}{$pack_num}; 
		  }
		 
	      }
	      else
	      {
		  ($not_verbose)or print "count_tot de $ind_tot $count_tot[$ind_tot] < $min_nb_of_waited_sigma_motifs\n";
	      }
	  }
	  
	  
	  
	  # if decallage to big, it consists in reinitialization
	  if( $decallage > $#count ){
	      $#count = -1;
	      $prev_subsp = '';
	      $read_keys_count--;
	      redo COUNTLOOP;
	  }
	  else
	  {
	      # decount of erases boxes
	      for my $decount( @count[0..$decallage - 1] )
	      { 
		  $count_summ -= $decount; 
		  if( $count_summ < 0 ){ die "$prog_tag [Error] countsumm negatif, line ".__LINE__."\n";}
	      }
	      
	      # real shifts of table boxes
	      @count[0..$#count - $decallage]     = @count[$decallage..$#count];
	      @count_tot[0..$#count - $decallage] = @count_tot[$decallage..$#count];
	      @count_tot[$#count - $decallage + 1..$#count] = (0)x $decallage;
	      
	      # new boxes to fill for counting
	      for my $i( $#count - $decallage + 1..$#count )
	      {
		  my $count_tmp = 0;
		  my $pack_subsp_v = pack('U',$subsp + $i - $indice_middle_spacershift_table);
		  
		  &COUNT_NB_SEQ($r_H_to_count_seq_by_spacer, \$count_tmp, $pack_subsp_v, $subsp + $i - $indice_middle_spacershift_table, $bool_IDtable);

		  $count[$i] = $count_tmp;

		  # print "on ajoute $count_tmp a count_summ qui valait $count_summ\n";
		  $count_summ += $count_tmp;
	      }
	  }
	  $prev_subsp = $subsp;
      }
      for my $ind(0..$#count_tot)
      {
	  if($count_summ > $count_tot[$ind]){ $count_tot[$ind] = $count_summ; }
      }
      
      # if we read the last key for spacement
      if( $read_keys_count == $tot_nb_keys )
      { 
	  for my $i_count(0..$#count)
	  {
	      if( ($count[$i_count] > 0)and($count_tot[$i_count] >= $min_nb_of_waited_sigma_motifs) )
	      {
		  my $pack_record_sp = pack('U',$subsp  - ( $indice_middle_spacershift_table - $i_count ));
		  if($bool_IDtable)
		  {
		      $r_valid_int_sp->{ $pack_record_sp } = 1;
		      # $r_valid_int_sp->{ $prev_subsp - ($indice_middle_spacershift_table - $i_count) } = 1;
		  }
		  else
		  {
		      ($not_verbose)or print "count_tot de $i_count $count_tot[$i_count] > $min_nb_of_waited_sigma_motifs\n";
		      &TREATMENT_COUNT_NB_SEQ($r_H_to_count_seq_by_spacer, $count_tot[$i_count], $num_subw1, $num_subw2, $pack_record_sp);
		  }
	      }
	      else
	      {
		  ($not_verbose)or print "count_tot de $i_count $count_tot[$i_count] < $min_nb_of_waited_sigma_motifs\n";
	      }
	  }
      }
      # else
      # {
	  # the spacement we treat is not the last, we do not have to record for the moment
      # }
      # print "FIN TRAITEMENT $table_trinuc[$num_subw1] $table_trinuc[$num_subw2] $subsp\n";
  }
} # end of sub TEST_NBSEQ_PER_INTERV
# ****************************************************************


# ?? INUTILE?
# sub REV_LENGTH_TRI_H(\@\@)
# {
  #  my ($r_input_array) = @_;
  #  @$r_input_array = map { $_->[1] }
  #  reverse sort { $a->[0] <=> $b->[0] }
  #  map { [ length $_, $_ ] } @$r_input_array;
  #  return $r_input_array;
# }
# &LENGTH_TRI_H(\@toto);

# **********************************

sub REMPLACE_W_PAR_ANTISLASHW(\$)
    {
	# word where we must replace 'n' by '\w'
    my ($r_re) = @_;
    # idem for 'n' and not 'w'
    my @w_pos = ();
    my $tmp_re = $$r_re;
    
    # necessary orelse we affect undefined value to $$r_re at the end
    ($$r_re =~ 'w')or return;

    while($$r_re =~ /(w+)/g)
    {
	# we record   pos of w+ match and length of w series
	push @w_pos, [pos($$r_re)-length($1), length($1)];
    }
    for my $iw(reverse 0..$#w_pos)
    {
	my $wregexp = '(\w{'.$w_pos[$iw][1].'})';
	($w_pos[$iw][1] == 1)and $wregexp = '(\w)';
	
	# print "partie precedente substr($tmp_re,0,$w_pos[$iw][0]) ";
	# print 'vaut ',substr($tmp_re,0,$w_pos[$iw][0]);
	# print "partie suivante substr($tmp_re, $w_pos[$iw][0]+$w_pos[$iw][1], length($tmp_re)-($w_pos[$iw][0]+ $w_pos[$iw][1])) ";
	# print 'vaut ',substr($tmp_re, $w_pos[$iw][0]+$w_pos[$iw][1], length($tmp_re)-($w_pos[$iw][0]+ $w_pos[$iw][1]))."\n";
	# print "on remplace $re par ";

	# on remplace les series de x 'w' par (\w{x}) et les 'w' simples par (\w)
	$tmp_re = substr($tmp_re,0,$w_pos[$iw][0]).$wregexp.substr($tmp_re, $w_pos[$iw][0]+$w_pos[$iw][1], length($tmp_re)-($w_pos[$iw][0]+ $w_pos[$iw][1]));
	# print "$tmp_re\n";
    }
    
    $$r_re = $tmp_re;
} # end of sub REMPLACE_W_PAR_ANTISLASHW(


# **********************************
sub REMPLACE_N_PAR_ANTISLASHW(\$)
{
    # word where we must replace 'n' by '\w'
    my ($r_re) = @_;
    # idem for 'n' and not 'w'
    my @w_pos = ();
    my $tmp_re = $$r_re;
    
    # necessary orelse we affect undefined value to $$r_re at the end
    ($$r_re =~ 'n')or return;

    while($$r_re =~ /(n+)/g)
    {
	# we record   pos of w+ match and length of w series
	push @w_pos, [pos($$r_re)-length($1), length($1)];
    }
    for my $iw(reverse 0..$#w_pos)
    {
	my $wregexp = '\w{'.$w_pos[$iw][1].'}';
	($w_pos[$iw][1] == 1)and $wregexp = '\w';
	
	# on remplace les series de x 'n' par \w{x} et les 'n' simples par \w
	$tmp_re = substr($tmp_re,0,$w_pos[$iw][0]).$wregexp.substr($tmp_re, $w_pos[$iw][0]+$w_pos[$iw][1], length($tmp_re)-($w_pos[$iw][0]+ $w_pos[$iw][1]));
    }
    
    $$r_re = $tmp_re;
} # end of sub REMPLACE_N_PAR_ANTISLASHW(
# **********************************

# **********************************
sub REMPLACE_N_PAR_POINT(\$)
{
    # word where we must replace 'n' by '\w'
    my ($r_re) = @_;
    # idem for 'n' and not '.'
    my @w_pos = ();
    my $tmp_re = $$r_re;
    
    # necessary orelse we affect undefined value to $$r_re at the end
    ($$r_re =~ 'n')or return;

    while($$r_re =~ /(n+)/g)
    {
	# we record   pos of w+ match and length of w series
	push @w_pos, [pos($$r_re)-length($1), length($1)];
    }
    for my $iw(reverse 0..$#w_pos)
    {
	my $wregexp = '.{'.$w_pos[$iw][1].'}';
	($w_pos[$iw][1] == 1)and $wregexp = '.';
	
	# on remplace les series de x 'n' par \w{x} et les 'n' simples par .
	$tmp_re = substr($tmp_re,0,$w_pos[$iw][0]).$wregexp.substr($tmp_re, $w_pos[$iw][0]+$w_pos[$iw][1], length($tmp_re)-($w_pos[$iw][0]+ $w_pos[$iw][1]));
    }
    
    $$r_re = $tmp_re;
} # end of sub REMPLACE_N_PAR_POINT(\$)
# **********************************


# **********************************FIN TRI RAPIDE LETTRE*********************************

# **** to know CPU and memory usage of the programme *************************************
# system('rm -rf MEM.txt');
# $_[0] is string to display
sub CAPTURE_COMPUTER_STAT($)
{
    # open(MEM,"> MEM.txt")or die "$!\n";
    open(TOP,"ps aux |" )or die "$prog_tag [Error] $!, line ".__LINE__."\n";
    # print MEM $_[0]."\n";
    print "******************\n".$_[0]."\n";
    # print MEM <TOP>;
    while( my $l = <TOP>)
    {
	($l =~ /(?:^USER|perl)/)and print $l;
    }
    close TOP;
    print "*******************************\n";
    # close MEM;
} # end of sub CAPTURE_COMPUTER_STAT($)
# ****************************************************************************************

# ***** subprog to know size of multihash or multilists
# sub SIZE_H_TABLE(\$)
# {
  #  my ($obj) = @_; 
  #  my $tmp_cpt = 0; 
    
  #  if(ref $$obj eq "SCALAR")
  #  {
	# return 1;
   # }
   # elsif(ref $$obj eq "REF")
   # {
	# return &SIZE_H_TABLE($$obj);
   # }
   # elsif(ref $$obj eq "HASH")
   # {
	# my $bool_final_H = 0;
	# while(my ($key,$val) = each(%$obj) )
	# {
	   # if(ref($val) eq "SCALAR"){ $bool_final_H = 1; }
	   # last; 
	# }
	# if($bool_final_H)
	# {
	  #  return scalar(keys %$obj);
	# }
	# else
	# {
	  #  while(my ($key,$val) = each(%$obj) )
	  #  {
		# $tmp_cpt += &SIZE_H_TABLE($val);
	   # }
	   # return $tmp_cpt;
	# }
    # }
    # elsif(ref $$obj eq "ARRAY")
    # {
	# if(ref($$obj[0]) eq "SCALAR")
	# { 
	  #  return ($#$obj +1);
	# }
	# else
	# {
	  #  for my $case(@$obj)
	  #  {
		# $tmp_cpt += &SIZE_H_TABLE($case);
	   # }
	   # return $tmp_cpt;
	# }
    # }
    
# }
# **************************************************************************************


sub RECORD_FINAL_RESULT($$\@$$\$$\@\$)
{
    my ($last_ind_seq_set_beg, $last_ind_seq_set_end, $r_tempo_spacer_bounds_ind_l_this_pos, $rapport_score, $rapport_de_vraisemblance, $r_regexp_motif, $nb_in_promot, $r_tmp_display , $r_string_let_pos_p) = @_;
     if(! $not_verbose)
     {
	print ' ' x length($$r_string_let_pos_p);
	print "SUB RECORD_FINAL_RESULT @_\n";
	print "rapport $rapport_score $$r_regexp_motif INTERV No ".(1+$#interv_int)." rapport_de_vraisemblance $rapport_de_vraisemblance\n";
	print "on ENREGISTRE last_ind_seq_set dans interv_int d'indice ".($#interv_int+1)."\n";
	
	print "SEQ IMPLIQUEES FINAL RESULT ln ".__LINE__."\n";
	for my $id_seq( $last_ind_seq_set_beg .. $last_ind_seq_set_end )
	{
	    my $real_ind_seq = $sorted_ind[$id_seq];
	    
	    
	    print "SEQ IMPLIQUEE $seq_sp[$sp][$real_ind_seq] :real_ind_seq $real_ind_seq, id_seq $id_seq\n"; # "dol1 $length length $length tempo_spacer_bounds @tempo_spacer_bounds\n";
	}
     }

    # TRI seq OK jusqu'ici, tri fait avant en fonction de l'espacement

    print "RECORD_FINAL_RESULT: rapport $rapport_score $$r_regexp_motif INTERV No ".(1+$#interv_int)." rapport_de_vraisemblance $rapport_de_vraisemblance\n";

   # print "motif enregistre $$r_regexp_motif\n";
    
    my $str_rapport_score            = sprintf "%2.2f ", $rapport_score;
    my $str_rapport_de_vraisemblance = sprintf "%2.2f ", $rapport_de_vraisemblance;

    # for french, german systems, to get english standards
	$str_rapport_score            =~ s/,/./g;
    $str_rapport_de_vraisemblance =~ s/,/./g;
    
    # $$r_string_let_pos_p             =~ s/_/ /g;

    # ($not_verbose)or print "SUB SEQ_SAME_SP @_\n";
    # ********************************************************************************************************
    # to align common motifs, we need to know the larger spacer to add gaps in other sequences
    my $longer_spacer = 0;
    my $length_same_begin = $length_same_begin_init;
    
    if($intitules_ind_prem_l !~ /^$deb_prev_intit[$sp]/ )
    {
		print {$fic_sp} "\n$intitules_ind_prem_l\n";
		$deb_prev_intit[$sp] = substr($intitules_ind_prem_l, 0, 30);
    } 
    
    print {$fic_sp} "\n";
    
    # print {$fic_sp} "MOTIF $$r_regexp_motif, R: $str_rapport_score (>= $rapports_display[$sp]), LRT is $str_rapport_de_vraisemblance (>= $POISSON_THRESHOLD, a = $alpha_likelyhood_ratio) nb seq ".($last_ind_seq_set_end - $last_ind_seq_set_beg +1).", $nb_in_promot in promot\n"; # , letter proba: $$r_string_let_pos_p\n";
	print {$fic_sp} join(' ',	"MOTIF $$r_regexp_motif, R:",
								$str_rapport_score,
								'(>=', $rapports_display[$sp],
								'), LRT is',
								$str_rapport_de_vraisemblance,
								'(>=', $POISSON_THRESHOLD,', a =',
								$alpha_likelyhood_ratio,
								') nb seq '.($last_ind_seq_set_end - $last_ind_seq_set_beg +1).", $nb_in_promot in promot\n"); # , letter proba: $$r_string_let_pos_p\n";
	
    
    # for all the sequences
    for my $l_not_sorted($last_ind_seq_set_beg.. $last_ind_seq_set_end)
    {
		# we get the first sorted indice to know which seq we have to display
		my $l = $sorted_ind[$l_not_sorted];
		$not_verbose or print "we record in file seq with real ind $l_not_sorted $seq_sp[$sp][$l], which would correspond to $l_not_sorted relative ind ln ".__LINE__."\n";
	
		print {$fic_sp} "$seq_sp[$sp][$l] $ID_f[$sp][$l]\n";
    }
    
    if(! $not_verbose)
    {
		print ' ' x length($$r_string_let_pos_p);
		print "FIN SUB RECORD_FINAL_RESULT @_\n";
    }

   # if($$r_regexp_motif eq "[at]aacc\\w{17,19}ccg"){
   #   print "BAD motif with not related seq line ".__LINE__."\n";
   #   exit;
   # }

} # end of sub RECORD_FINAL_RESULT($$\@$$\$$\@\$)


# ***************************************************************************
# DEAL WITH EVALUATION AND RUN OF RECURSIVE PROCESS
# ***************************************************************************
# return 0 if no interesting motif is found in larger set of seq, 1 otherelse (used to know if
# we have to go on extending
sub EVALUATION($$$\@$\@\@\$$$\@\@$$\@){
    
    my ($last_ind_seq_set_beg_l, $last_ind_seq_set_end_l, $bool_split_in_record_not_done, $r_tempo_spacer_bounds_ind_l_this_pos, $pos_proba_min_0, $r_proba_ext, $r_border_words, $string_let_pos_p, $bool_not_recorded_interval_l, $r_regexp_motif, $r_low_proba, $r_global_cpt_score_pos_proba_min, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice, $r_last_ind_seq_set) = @_;

    # print "sub EVALUATION at ",scalar(localtime(time())),"\n";

    if(! $not_verbose){
      print "SUB Evaluation @_\n";
      print ' ' x length($string_let_pos_p);
      print "SUB Evaluation @_\n";
    }

    my $returned_val; # returned val by EVALUATION (0 interesting motif found, 0 N
    my @tmp_display = ();
    
    # print "last_ind_seq_set_beg_l : $last_ind_seq_set_beg_l, ";
    # print "last_ind_seq_set_end_l : $last_ind_seq_set_end_l\n";
    
    # *******************************  SELECTION CRITERIA  *********************************************
    my $rapport_score      = -1;
    my $nb_in_promot       = 0;
    my $nb_in_whole_2_sens = 0;
    
    # if extension of motif is to long, this is probably not a TFBS and statistics are false because of length of the motif, we do not need to compute
    # DECLARER filtre_nb_l_added_matrice_mini EN PARAMETRE ?? !!!!!!!!!!!
    my $nbl = $nb_l_added_matrice + $initial_total_nb_boxes_letters; 
    if(($nbl >= $filtre_nb_l_added_matrice_mini)and($nbl <= $filtre_nb_l_added_matrice_maxi)){
	&RAPPORT_CRITERIA( $r_regexp_motif, \$rapport_score, \$nb_in_promot, \$nb_in_whole_2_sens);
	($not_verbose)or print "rapport $rapport_score $$r_regexp_motif INTERV No ".(1+$#interv_int)." for r_regexp_motif $$r_regexp_motif\n";
    }
    else{

	($not_verbose)or print "EXTENSION TO short FOR $$r_regexp_motif , $nb_in_promot in promot\nnb_l_added_matrice + initial_total_nb_boxes_letters ".($nb_l_added_matrice + $initial_total_nb_boxes_letters)." < filtre_nb_l_added_matrice_mini $filtre_nb_l_added_matrice_mini\n";
	return 0;
    }
    
    # tmp array to save sorted_ind array because every recorded data are given for sorted_ind with the order found here (other sorts can change the order of seq)
    # if no result is found for data after other sorts, we will need to recover sorted_ind
    my @memory_sort         = ();
    my @memory_pos_tri_tmp  = ();

    if(($rapport_score >= $rapports[$sp])and($bool_not_recorded_interval_l)){

	my ($P_whole, $P_upstr, $max_motif_length) = (0)x2;

	&max_motif_length_and_P_MM0($r_regexp_motif, \$max_motif_length, 1, \$P_whole, \$P_upstr);
	# die "P_whole, P_upstr: $P_whole, $P_upstr\n";	

	my $theor_length_of_upstr_seq = $upstr_nt[$sp] - $nb_upstr_seq[$sp] * $max_motif_length;
	my $theor_length_of_genome_seq = ($total_nt[$sp] - $nb_whole_seq[$sp] * $max_motif_length)  << 1;
	my $nb_total = $nb_in_promot + $nb_in_whole_2_sens;
	my $pi = ($theor_length_of_upstr_seq * $P_upstr) / ( $theor_length_of_upstr_seq * $P_upstr + $theor_length_of_genome_seq * $P_whole);

	# ($not_verbose)or print "For R and LRT computing, nb_in_promot $nb_in_promot nb_in_whole_2_sens $nb_in_whole_2_sens\n";

	# *********************************** OLD **********************************************
	# my $N_sur_l = ($nb_in_whole_2_sens + $nb_in_promot)/($theor_length_of_upstr_seq + $theor_length_of_genome_seq);

	# ($not_verbose)or print "first log ".(log(($nb_in_promot/$theor_length_of_upstr_seq)/($N_sur_l)))." second log ".(log(($nb_in_whole_2_sens/$theor_length_of_genome_seq)/($N_sur_l)))."\n";
	# my $rapport_de_vraisemblance = ($nb_in_promot * log(($nb_in_promot/$theor_length_of_upstr_seq)/($N_sur_l))
	#				+ $nb_in_whole_2_sens * log(($nb_in_whole_2_sens/$theor_length_of_genome_seq)/($N_sur_l))
	#				) * 2;
	# *********************************** OLD **********************************************
	# *********************************** NEW MM0 ******************************************
	my $rapport_de_vraisemblance = ($nb_in_promot * log(($nb_in_promot/$nb_total)/$pi)
					+ $nb_in_whole_2_sens * log(($nb_in_whole_2_sens/$nb_total)/(1-$pi))
					) * 2;



	# *********************************** NEW MM0 ******************************************

	# print "theor_length_of_upstr_seq $theor_length_of_upstr_seq, theor_length_of_genome_seq $theor_length_of_genome_seq, N_sur_l $N_sur_l\n";

	# ($not_verbose)or print "rapport_de_vraisemblance $rapport_de_vraisemblance is:\n($nb_in_promot * log(($nb_in_promot/$theor_length_of_upstr_seq)/($N_sur_l)) + $nb_in_whole_2_sens * log(($nb_in_whole_2_sens/$theor_length_of_genome_seq)/($N_sur_l))) * 2\n";
	# *******************************  END SELECTION CRITERIA  ***************************************** #
	
	
	
	if($rapport_de_vraisemblance >= $POISSON_THRESHOLD){
	    
	    
	  # ****************** #
	  
	  # WHEN we have an interesting motif, we have to verify if it is not composed of two, with really distinct spacers, like:
	  # ggaat\w{10, 23}gtt could be composed of ggaat\w{10,12}gtt and ggaat\w{20,23}gtt... it is very different

	  # If it is the case, we evaluate each set of seq and record them if related motif is interesting (run of &evaluate),,
	  # otherelse, we record global sequence set found as interesting

	    # to apply filter only in first call
	    if($bool_split_in_record_not_done){

		$bool_split_in_record_not_done = 0;
		# ****************** we delete seq corresponding to orphelin spacer and eventually
		# hash to count seq by spacer (with spacershift allowed) to split seq set in subset
		my %H_to_count_seq_by_spacer = ();
		
		# ($not_verbose)or print "ind_beg_to_sort $ind_beg_to_sort ind_end_to_sort $ind_end_to_sort\n";
		
		# count of number of seq foreach sp (sort, 'delete', split seq set)

		# record of usefull sorted_ind part for eventually recover it later
		push @memory_sort, @sorted_ind[$last_ind_seq_set_beg_l .. $last_ind_seq_set_end_l];
		@memory_pos_tri_tmp = @{ &clone($pos_tri_tmp[$sp]) };
		$not_verbose or print "clone pos_tri_tmp line ".__LINE__."\n";

		($#{ $pos_tri_tmp[$sp] } != $#memory_pos_tri_tmp)and die "$prog_tag [Error] clone does not work line ".__LINE__."\n";
		# we store every real ind of seq in hash table H_to_count_seq_by_spacer for a given spacer (key)
		for my $id_seq( $last_ind_seq_set_beg_l .. $last_ind_seq_set_end_l ){

		    my $real_ind_seq = $sorted_ind[$id_seq];
		    
		    if(! $not_verbose){

			print "SEQ IMPLIQUEE $seq_sp[$sp][$real_ind_seq] d'ind $real_ind_seq\n"; # "dol1 $length length $length tempo_spacer_bounds @tempo_spacer_bounds\n";
			# print "r_pos_tri de 0 ".$pos_tri[$sp][$real_ind_seq][0]."\n";
			# print "r_pos_tri de 1 ".$pos_tri[$sp][$real_ind_seq][1]."\n";
		    }
		    
		    my $spacer_l = pack('U', $pos_tri[$sp][$real_ind_seq][1] - $pos_tri[$sp][$real_ind_seq][0] - $LOCAL_WIDTH_OF_L_HBOXES);

		    # HBOXES EN COURS
		    # push @{ $H_to_count_seq_by_spacer{$spacer_l} }, $id_seq;
		    push @{ $H_to_count_seq_by_spacer{$spacer_l} }, $real_ind_seq;
		    # ($not_verbose)or print "on a mis $id_seq d'espace ".unpack('U',$spacer_l)." dans @{ $H_to_count_seq_by_spacer{$spacer_l} }\n";
		}
		
		# OK
		if(! $not_verbose){

		    print "*************\nverif H_to_count_seq_by_spacer H\n";
		    for my $key(keys %H_to_count_seq_by_spacer)
		    {
			print "spacer ".unpack('U',$key).", seq indices: @{$H_to_count_seq_by_spacer{$key}}\n";
		    }
		    print "END verif H_to_count_seq_by_spacer H\n*************\n";
		    
		}
		
		my $indice_middle_spacershift_table = int($#spacer_shift_fic / 2); # subscript of the middle of the array
		# 3 -> 1
		# 5 -> 2
		my %valid_int_sp = (); 
		
		# count number of distinct seq involved for a given spacer (by taking into account spacershift)
		&TEST_NBSEQ_PER_INTERV(\%H_to_count_seq_by_spacer, $indice_middle_spacershift_table, \%valid_int_sp); # TEST valid_int_sp
		# if(!$not_verbose)
		# {
		  #  if(scalar(keys %valid_int_sp) == 0)
		  #  {
			# print "valid_int_sp H do not exists!!\n";
		  #  }
		  #  else
		  #  {
			# print "valid_int_sp\n";
			# while(my ($clef,$val) = each %valid_int_sp)
			# {
			  #  print "subw1 $border_words[0] subw2 $border_words[1]: clef (space) ".unpack('U',$clef)."\n";
			# }
			# print "FIN valid_int_sp\n";
			
			
		    # }
		# }
		
		# BEGINNING OF SPLIT***********************************************************
		# split into separate interv, 'delete' of seq with ininteresting spacer
		my @sort_ind_to_sort_l   = (); # sorted sorted_ind array (interesting subscripts)
		my @sort_ind_to_sort_bad = (); # sorted sorted_ind array (uninteresting subscripts)
		my @ind_to_split         = (); # store last subscript of interesting seq set
		
		my $previous_interv_valid = 0;
		
		# for every met spacer sorted by increasing order...
		for my $sorted_sp(sort unpack_numerique_U keys %H_to_count_seq_by_spacer){

		    my $unpack_sorted_sp = unpack('U', $sorted_sp);
		    ($not_verbose)or print "sorted_sp $unpack_sorted_sp\n";

		    # if we have at less a second interesting interval
		    
		    if(exists $valid_int_sp{$sorted_sp}){

			($not_verbose)or print "valid space $unpack_sorted_sp\n";
			
			# if we have an interesting seq set (enough seq), we store them into sort_ind_to_sort_l (so this is sorted)
			push @sort_ind_to_sort_l, @{ $H_to_count_seq_by_spacer{$sorted_sp} };

			# print "on met @{ $H_to_count_seq_by_spacer{$sorted_sp} } ds sort_ind_to_sort_l pour unpack_sorted_sp\n";
			if($previous_interv_valid){
			    
			  # if previous met sequences set was recorded, this means that previous and current one
			  # have to be merged, so we replace previous recorded last subscript by a new one corresponding
			  # to those of last seq of current sequences set
			  $ind_to_split[$#ind_to_split] = $last_ind_seq_set_beg_l + $#sort_ind_to_sort_l;
			}
			else{

			    # we record because it is at less the second interesting intervalle and we have to run recursiverly the programme on... 
			    push @ind_to_split, $last_ind_seq_set_beg_l + $#sort_ind_to_sort_l;
			    ($not_verbose)or print "limit of the indice to sort, recorded in ind_to_split $ind_to_split[$#ind_to_split]\n";
			}
			# we record last space used
			# print "valid_int_sp of $unpack_sorted_sp ".$valid_int_sp{$sorted_sp}."\nmust not be empty\n"; # ??
			$previous_interv_valid = 1;
		    }
		    else{

			if(! $not_verbose){
			    print "nb ds tab: ".$#{$H_to_count_seq_by_spacer{$sorted_sp}}." NOT INTERESTING for spacer ".unpack('U',$sorted_sp)." \n";
			    print "ref: $H_to_count_seq_by_spacer{$sorted_sp}\ntab: ".@{$H_to_count_seq_by_spacer{$sorted_sp}}."\n";
			}
			
			# print "on met @{ $H_to_count_seq_by_spacer{$sorted_sp} } ds sort_ind_to_sort_bad pour unpack_sorted_sp\n";
			push @sort_ind_to_sort_bad, @{ $H_to_count_seq_by_spacer{$sorted_sp} };
			$previous_interv_valid = 0;
		    }
		}
		# we have to consider sequences subscripts from 0 to last_ind_seq_set_beg_l-1 as not interesting because they are not treated
		# normally by the run of this sub
		push @sort_ind_to_sort_bad, @sorted_ind[0..$last_ind_seq_set_beg_l-1];
		# print "on met @sorted_ind[0..$last_ind_seq_set_beg_l-1] d'indices 0..".($last_ind_seq_set_beg_l-1)." ds sort_ind_to_sort_bad\n";

		# idem for sequences subscripts after last treated seq
		push @sort_ind_to_sort_bad, @sorted_ind[$last_ind_seq_set_end_l+1..$#sorted_ind];
		# print "on met @sorted_ind[$last_ind_seq_set_end_l+1..$#sorted_ind]  d'indices ".($last_ind_seq_set_end_l+1)." $#sorted_ind ds sort_ind_to_sort_bad\n";

		($not_verbose)or print "sorted_ind @sorted_ind AV REMPLACE PAR\n";

		# print "sorted_ind 0..$#sort_ind_to_sort_l recoit @sort_ind_to_sort_l\n";
		@sorted_ind[0..$#sort_ind_to_sort_l] = @sort_ind_to_sort_l;
		# print "sorted_ind ".($#sort_ind_to_sort_l+1)."..$#sorted_ind recoit @sort_ind_to_sort_bad\n";
		@sorted_ind[$#sort_ind_to_sort_l+1..$#sorted_ind] = @sort_ind_to_sort_bad;

		($not_verbose)or print "sorted_ind @sorted_ind AP\n";
		
		# if we have only 1 value, it corresponds to the first intervalle treated after, we do not have other intervals to treat,
		# so we do not need to enter hear (otherelse we only use indice of first intervalle to know from which indice we have to 
		# consider the second intervale (the first we have to treat hear in fact)
		if($#ind_to_split > 0){

		    ($not_verbose)or print "****************************************\non scinde l'intervalle en plusieurs intervalles\n";
		    
		    my $nb_l_added_matrice_l = $initial_total_nb_boxes_letters + $nb_l_added_matrice - 1;

		    my $bool_interesting_sub_res = 0; # bool to know if we have found an interesting subset of sequences with lower spacer variation

		    # treatment of intervalles of sequences corresponding to different spacer ranges for interesting sequences
		    for my $i_to_sort(0..$#ind_to_split-1){
				    
		      my @spacer_bound_splited_interv = (100,0);

		      ($not_verbose)or  print "******************************\nverif SEQ for $i_to_sort interval in SUB SORT, seq No $ind_to_split[$i_to_sort] from ".($ind_to_split[$i_to_sort]+1)." to ".($ind_to_split[$i_to_sort+1])."\n";
		      
		      # we record in spacer_bound_splited_interv lowest and highest spacer between our 2 'words'
		      for my $i($ind_to_split[$i_to_sort]+1..$ind_to_split[$i_to_sort+1]){
			
			my $r_i = $sorted_ind[ $i ];
			# print $seq_sp[$sp][$r_i]." id $r_i\n";  
			
			
			# pos at right of -35 box #**_____*** 
			# $pos_tri[$sp][$exchanged_real_ind][0]; 
			# pos at left of -10 box ***_____#**
			# $pos_tri[$sp][$exchanged_real_ind][1]; 
			my $spacer_split_interv = $pos_tri[$sp][$r_i][1] - $pos_tri[$sp][$r_i][0] - $LOCAL_WIDTH_OF_L_HBOXES;
			
			# lower spacer between seeds
			( $spacer_split_interv < $spacer_bound_splited_interv[0] )and
			  $spacer_bound_splited_interv[0] = $spacer_split_interv;
			
			# higher spacer between seeds
			(  $spacer_split_interv > $spacer_bound_splited_interv[1] )and
			  $spacer_bound_splited_interv[1] = $spacer_split_interv;
			
		      }
		      ($not_verbose)or print "\nverif seq for $i_to_sort interval in SUB SORT\n******************************\n";
		      ($ind_to_split[$i_to_sort+1] - $ind_to_split[$i_to_sort] < $nb_mini_occ)and die "$prog_tag [Error] PB: Not enough sequence in interesting interv line ".__LINE__."\n";
		      
		      # print "spacer_bound_splited_interv @spacer_bound_splited_interv\n";
		      
		      &INIT_PROBA_TAB($r_proba_ext, 
				      $r_border_words, 
				      $ind_to_split[$i_to_sort+1]- $ind_to_split[$i_to_sort]
			  );
		      
		      my $regexp_motif_l = '';
		      $not_verbose or print "MATRICE_CREATE2 line ".__LINE__." nb_seq ".($ind_to_split[$i_to_sort+1]- $ind_to_split[$i_to_sort])."\n";
		      &MATRICE_CREATE2($r_border_words, 
				       $r_low_proba->[$pos_proba_min_0], 
				       \$regexp_motif_l, 
				       \@spacer_bound_splited_interv, 
				       $r_global_cpt_score_pos_proba_min, 
				       $ind_to_split[$i_to_sort+1]- $ind_to_split[$i_to_sort], 
				       $nb_of_intern_l_added_to_boxes, 
				       \$nb_l_added_matrice_l, 
				       $pos_proba_min_0);
		      
		      # Here, we consider individually seq set with letter of low proba at pos_proba_mean, proba are not reinitialized as they do not change
		      # but evaluation is different because number of seq involved is different
		      
		      # print "regexp_motif_l $regexp_motif_l\n";
		      
		      if($nb_l_added_matrice_l != $nb_l_added_matrice){

			$not_verbose or print "On lance EVAL line ".__LINE__."\n";
			my @local_ind_to_split = @ind_to_split[$i_to_sort..$i_to_sort+1]; # limit next recursives runs to this sole interval (otherelse, would take into account other split intervals)
			$bool_interesting_sub_res = $bool_interesting_sub_res or &EVALUATION($ind_to_split[$i_to_sort]+1, $ind_to_split[$i_to_sort+1], 0, \@spacer_bound_splited_interv, $pos_proba_min_0, $r_proba_ext, $r_border_words, $string_let_pos_p, $bool_not_recorded_interval_l, \$regexp_motif_l, $r_low_proba, $r_global_cpt_score_pos_proba_min, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice_l, \@local_ind_to_split);
			
		      }
		      else{
			$not_verbose or print "No new letter in MATRICE_CREATE2, we return!!! line ".__LINE__."\n"; 
		      }
		    }

		    if($bool_interesting_sub_res){
		      $bool_not_recorded_interval_l = 0;
		      if($#memory_sort != -1){
			@sorted_ind[$last_ind_seq_set_beg_l, $last_ind_seq_set_end_l] = @memory_sort;
			$pos_tri_tmp[$sp] = &clone(\@memory_pos_tri_tmp);
			$not_verbose or print "clone pos_tri_tmp line ".__LINE__."\n";
			($#{ $pos_tri_tmp[$sp] } != $#memory_pos_tri_tmp)and die "$prog_tag [Error] clone does not work line ".__LINE__."\n";
		      }
		      # print "ICI: ind_to_split @ind_to_split\n";

		      $not_verbose or print "ON ENREGISTRE LE RESULTAT FINAL line ".__LINE__."\n";		
		      &RECORD_FINAL_RESULT($last_ind_seq_set_beg_l, $last_ind_seq_set_end_l, $r_tempo_spacer_bounds_ind_l_this_pos, $rapport_score, $rapport_de_vraisemblance, $r_regexp_motif, $nb_in_promot, \@tmp_display , \$string_let_pos_p);
		      $returned_val = 1;
		    }

		    # treatment of last intervalle of sequences corresponding to different spacer ranges for interesting sequences
		    # ($not_verbose)or print "on lance SORT_IND indices $ind_to_split[$#ind_to_split], $#sort_ind_to_sort_l DANS SPLIT\n";

		    if(! $not_verbose){

			print "******************************\nverif SEQ for last interval in SUB SORT, seq No $ind_to_split[$#ind_to_split] to $#sort_ind_to_sort_l\n";
			for my $i($ind_to_split[$#ind_to_split]..$#sort_ind_to_sort_l){

			    my $r_i = $sorted_ind[ $i ];
			    print $seq_sp[$sp][$r_i]." id $r_i\n";  
			}
			print "\nverif seq for last interval in SUB SORT\n******************************\n";
		    }
		    
		}
		elsif($last_ind_seq_set_end_l - $last_ind_seq_set_beg_l + 1 >= $nb_mini_occ ){

		    $bool_not_recorded_interval_l = 0;
		    if($#memory_sort != -1){
		      @sorted_ind[$last_ind_seq_set_beg_l, $last_ind_seq_set_end_l] = @memory_sort;
		      $pos_tri_tmp[$sp] = &clone(\@memory_pos_tri_tmp);
		      $not_verbose or print "clone pos_tri_tmp line ".__LINE__."\n";
		      ($#{ $pos_tri_tmp[$sp] } != $#memory_pos_tri_tmp)and die "$prog_tag [Error] clone does not work line ".__LINE__."\n";
		    }
		    # print "ICI: ind_to_split @ind_to_split\n";
		    $not_verbose or print "ON ENREGISTRE LE RESULTAT FINAL line ".__LINE__."\n";		
		    &RECORD_FINAL_RESULT($last_ind_seq_set_beg_l, $last_ind_seq_set_end_l, $r_tempo_spacer_bounds_ind_l_this_pos, $rapport_score, $rapport_de_vraisemblance, $r_regexp_motif, $nb_in_promot, \@tmp_display , \$string_let_pos_p);
		    $returned_val = 1;
		}
	    }
	    # END OF SPLIT***********************************************************

	    # split already done, nb seq sufficient, we record result
	    elsif($last_ind_seq_set_end_l - $last_ind_seq_set_beg_l + 1 >= $nb_mini_occ ){

		$bool_not_recorded_interval_l = 0;
		if($#memory_sort != -1){
		  @sorted_ind[$last_ind_seq_set_beg_l, $last_ind_seq_set_end_l] = @memory_sort;
		  $pos_tri_tmp[$sp] = &clone(\@memory_pos_tri_tmp);
		  $not_verbose or print "clone pos_tri_tmp line ".__LINE__."\n";
		  ($#{ $pos_tri_tmp[$sp] } != $#memory_pos_tri_tmp)and die "$prog_tag [Error] clone does not work line ".__LINE__."\n";
		}
		$not_verbose or print "ON ENREGISTRE LE RESULTAT FINAL line ".__LINE__."\n";		
		&RECORD_FINAL_RESULT($last_ind_seq_set_beg_l, $last_ind_seq_set_end_l, $r_tempo_spacer_bounds_ind_l_this_pos, $rapport_score, $rapport_de_vraisemblance, $r_regexp_motif, $nb_in_promot, \@tmp_display, \$string_let_pos_p );
		$returned_val = 1;
	    }
	}
	else{

	    ($not_verbose)or print "rapport $rapport_score high but rapport_de_vraisemblance $rapport_de_vraisemblance low (<$POISSON_THRESHOLD) for $$r_regexp_motif\n";
	    $bool_not_recorded_interval_l = 1;
	    $returned_val = 0;
	}
    }
    elsif($bool_not_recorded_interval_l){

	# rapport de vraisemblance > poisson threshold, we have to extend
	($not_verbose)or print  "rapport $rapport_score < $rapports_display[$sp]!!!!\n";

	# NOOOOOOO: memory_sort has not been recorded if we are here because no other sort of sorted_ind was done....useless
	# ($#memory_sort == -1)or @sorted_ind[$last_ind_seq_set_beg_l, $last_ind_seq_set_end_l] = @memory_sort;
	# $bool_not_recorded_interval_l = 1;
	$returned_val = 0;
    }
    
    # conditions to stop searches:
    #  interval not marked as interesting and nb of involved seq inf to minimal number of involved
    #  interesting interval but extension involves to few seq
    # print "last_ind_seq_set_end_l $last_ind_seq_set_end_l,";
    # print " last_ind_seq_set_beg_l $last_ind_seq_set_beg_l\n";
    if( ($bool_not_recorded_interval_l)and($last_ind_seq_set_end_l - $last_ind_seq_set_beg_l +1  >= $nb_mini_occ)){
    # if( ($last_ind_seq_set_end_l - $last_ind_seq_set_beg_l +1  >= $nb_mini_occ)or
	# ((not $bool_not_recorded_interval_l) and ($last_ind_seq_set_end_l - $last_ind_seq_set_beg_l +1  >= $nb_mini_of_seq_by_interv))	){

	my @bool = (0,0,0,0);   
	$bool[ $pos_proba_min_0 ] = 1;	
	
	# we copy current to avoid effects on other sort prog launched 

	# usefull???
	my @proba_ext = ();
	# @proba_ext = @{ &clone( \@$r_proba_ext ) }; # $proba_ext[ pos ][ lettre(numero) ][ k ]: k de Cnk	
	
	# &INIT_PROBA_TAB($r_proba_ext, \@multi_border_words, $last_ind_seq_set[$ind_l_this_pos+1] - $last_ind_seq_set[$ind_l_this_pos] -1, \@bool);
	
	# do nothing if called from treatment_one_low_proba (border_words already transform but we have extended of one letter, which is not the case when

    # FOR GOING ON extensions in different way &eval for each low proba and &recup to extend *************************************
	my @extended_border_words = @$r_border_words;
	my $lettre = '';
	
	if($bool_degenerescence_max_motif){
	  # TO TREAT ******************************************************* ??

	  die "$prog_tag [Error] not treated; we can not be here line ".__LINE__."\n";
	    # NEEDS MODIFICATIONS FOR PROBA EXTENSION TO WORK: MM3 on [gt]acc impossible for example for
	    # the moment because we need fixed letter as gacc (we extract trinuc key thanks substr(border,0,3)
	    # for instance ??
	    # **************************************************************************************************
	    for(@{ $r_low_proba->[$pos_proba_min_0] }){
		$lettre .= $rev_H_l[ $r_low_proba->[$pos_proba_min_0][$_][0] ]; # sorted tab, so letter with the lowest proba in first pos
		# print "r_low_proba_pos_proba_min $r_low_proba->[$pos_proba_min_0][$_][0] doit correspondre a la lettre $lettre ln ".__LINE__."\n";	    
		
	    }

	    (length($lettre) > 1)and $lettre = '['.$lettre.']'; 

	    # CAUTION nb_of_intern_l_added_to_boxes increased for futur calls &eval and &recup... DO NOT ADD other sort or treatments
	    if   ( $pos_proba_min_0 == 0 ){ $extended_border_words[0] = $lettre.$r_border_words->[ 0 ]; }
	    elsif( $pos_proba_min_0 == 2 ){ $extended_border_words[1] = $lettre.$r_border_words->[ 1 ]; $nb_of_intern_l_added_to_boxes++; }
	    # extension on the right
	    elsif( $pos_proba_min_0 == 1 ){ $extended_border_words[0] = $r_border_words->[ 0 ].$lettre; $nb_of_intern_l_added_to_boxes++; }
	    elsif( $pos_proba_min_0 == 3 ){ $extended_border_words[1] = $r_border_words->[ 1 ].$lettre; }

	    # PB: MM3 with boxes like [ac]gg to extend ??

	    ($not_verbose)or print "On lance RECUP_LETTRE_SORT @bool sur les sequences indices $last_ind_seq_set_beg_l $last_ind_seq_set_end_l, dans sub SORT avec border_words @extended_border_words bool_not_recorded_interval_l ".(($#interv_int == -1)or($last_ind_seq_set_beg_l > $interv_int[$#interv_int]))." line ".__LINE__."\n";
	    
	    # if we treat first interval (corresponding to the most overrepresented letter), we extend
	    # ($last_ind_seq_set_beg_l > $interv_int[$#interv_int]) means intervals has not been recorded
	    # &RECUP_LETTRE_SORT($last_ind_seq_set_beg_l, $last_ind_seq_set_end_l, @bool , (($#interv_int == -1)or($last_ind_seq_set_beg_l > $interv_int[$#interv_int])), \@extended_border_words_call, \@proba_ext, $bool_split_in_record_not_done);
	    
	    # + 1 is the letter which will be added thanks RECUP... useless as created regexp motif will be to long
	    if($nb_l_added_matrice + $initial_total_nb_boxes_letters + 1 >= $filtre_nb_l_added_matrice_maxi){
		
		($not_verbose)or print "EXTENSION TO LONG FOR futur regexp_motif, we do not need to get letters in RECUP , nb_l_added_matrice + initial_total_nb_boxes_letters ".($initial_total_nb_boxes_letters + $nb_l_added_matrice)." >= filtre_nb_l_added_matrice_maxi $filtre_nb_l_added_matrice_maxi\n";

		return 0;
	    }
	    
	    $not_verbose or print "lancement RECUP string_let_pos_p $string_let_pos_p, bool_not_recorded_interval_l $bool_not_recorded_interval_l, line ".__LINE__."\n";
	    
	    # $nb_l_added_matrice++;
	    $not_verbose or print "On lance RECUP $last_ind_seq_set_beg_l, $last_ind_seq_set_end_l, @bool , $bool_not_recorded_interval_l, \@extended_border_words, \@proba_ext, $bool_split_in_record_not_done, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice:line ".__LINE__."\n";
	    die "$prog_tag [Error] we must not be here: case not treated!!!!!!!, degenerated motif line ".__LINE__."\n";
	    &RECUP_LETTRE_SORT($last_ind_seq_set_beg_l, $last_ind_seq_set_end_l, @bool , $bool_not_recorded_interval_l, \@extended_border_words, \@proba_ext, $bool_split_in_record_not_done, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice);
	    # **************************************************************************************************	
	    # END TO TREAT ******************************************************* ??    
	}
	else{
	  # we need to get initial order to have coherent data
	  if($#memory_sort != -1){
	    if($last_ind_seq_set_end_l - $last_ind_seq_set_beg_l != $#memory_sort){ die "$prog_tag [Error] PB: we affect table with size which do not corresponds line ".__LINE__." (".($last_ind_seq_set_end_l - $last_ind_seq_set_beg_l)." diff de $#memory_sort)\n"; }
	    @sorted_ind[$last_ind_seq_set_beg_l, $last_ind_seq_set_end_l] = @memory_sort;
	    $pos_tri_tmp[$sp] = &clone(\@memory_pos_tri_tmp);
	    if(not $not_verbose){
	      print "clone pos_tri_tmp line ".__LINE__."\n";
	      print "we recover memory_sort @memory_sort line ".__LINE__."\n";
	    }
	    ($#{ $pos_tri_tmp[$sp] } != $#memory_pos_tri_tmp)and die "$prog_tag [Error] clone does not work line ".__LINE__."\n";
	  }

	  # print "r_last_ind_seq_set @$r_last_ind_seq_set\n";

	  $not_verbose or &verif_low_prob($r_low_proba, $pos_proba_min_0, __LINE__);
	   
	  # for my $ilp(0..$#{ $r_low_proba->[$pos_proba_min_0] }){
	  for my $ilp(0..$#$r_last_ind_seq_set-1){
	    # $#{ $r_low_proba->[$pos_proba_min_0] } or last;

	    # print "INTERVALLE ilp $ilp";
	    ($r_last_ind_seq_set->[$ilp+1]- $r_last_ind_seq_set->[$ilp] < $nb_mini_occ)and next;
	    # **************************************************************************************************	    
	    # evaluation of motif corresponding to only one letter at position of low_proba (and not every low_proba letters as previously)
	    
	    # print "*************************\nlow_proba  r_low_proba de $pos_proba_min_0 $ilp (ilp):  $r_low_proba->[$pos_proba_min_0][$ilp][0], $r_low_proba->[$pos_proba_min_0][$ilp][1] doit correspondre aux seq\n";
	    
	    my $bool_interess_res = 0;
	    
	    # if there is only one low_proba, we are sur motif has already be tested, so we don't test
	    
	    # if($#{ $r_low_proba->[$pos_proba_min_0] } > 0){		       	      
	    # we get range of spacer
	    my @spacer_bound_set_of_seq = (100,0);
	    
	    for my $i($r_last_ind_seq_set->[$ilp]+1..$r_last_ind_seq_set->[$ilp+1]){
	      
	      my $r_i = $sorted_ind[ $i ];
	      # print $seq_sp[$sp][$r_i]." id $r_i\n";  
	      
	      # pos at right of -35 box #**_____*** 
	      # $pos_tri[$sp][$exchanged_real_ind][0]; 
	      # pos at left of -10 box ***_____#**
	      # $pos_tri[$sp][$exchanged_real_ind][1]; 
	      my $spacer_split_interv = $pos_tri[$sp][$r_i][1] - $pos_tri[$sp][$r_i][0] - $LOCAL_WIDTH_OF_L_HBOXES;
		
	      # lower spacer between seeds
	      ( $spacer_split_interv < $spacer_bound_set_of_seq[0] )and
		  $spacer_bound_set_of_seq[0] = $spacer_split_interv;
	      
	      # higher spacer between seeds
	      (  $spacer_split_interv > $spacer_bound_set_of_seq[1] )and
		$spacer_bound_set_of_seq[1] = $spacer_split_interv;
	      
	    }
	    
	    #print "*************************\n";
	    
	    &INIT_PROBA_TAB($r_proba_ext, $r_border_words, $r_last_ind_seq_set->[$ilp+1]-$r_last_ind_seq_set->[$ilp]);
	    
	    my $regexp_motif_l = '';
	      
	    # my ($last_ind_seq_set_beg_l, $last_ind_seq_set_end_l, $bool_split_in_record_not_done, $r_tempo_spacer_bounds_ind_l_this_pos, $pos_proba_min_0, $r_proba_ext, $r_border_words, $string_let_pos_p, $bool_not_recorded_interval_l, $r_regexp_motif, $r_low_proba, $r_global_cpt_score_pos_proba_min, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice, $r_last_ind_seq_set)
	    
	    my $nb_l_added_matrice_l = $nb_l_added_matrice;
	    
	    $not_verbose or print "matrice create avec low_proba  r_low_proba de $pos_proba_min_0 $ilp (ilp):  $r_low_proba->[$pos_proba_min_0][$ilp][0], $r_low_proba->[$pos_proba_min_0][$ilp][1]\n";
	    $not_verbose or print "On lance MATRICE_CREATE2_1L_AT_1POS @$r_border_words, $r_low_proba->[$pos_proba_min_0][$ilp], $regexp_motif_l, @spacer_bound_set_of_seq, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice_l, $pos_proba_min_0\n";
	    &MATRICE_CREATE2_1L_AT_1POS($r_border_words, $r_low_proba->[$pos_proba_min_0][$ilp], \$regexp_motif_l, \@spacer_bound_set_of_seq, $nb_of_intern_l_added_to_boxes, \$nb_l_added_matrice_l, $pos_proba_min_0);
	    
	    if($nb_l_added_matrice_l != $nb_l_added_matrice){
	      $not_verbose or print "On lance EVAL line ".__LINE__."\n";
	      my @local_ind_to_split = @$r_last_ind_seq_set[$ilp..$ilp+1]; # limit next recursives runs to this sole interval (otherelse, would take into account other split intervals)
	      my @local_proba_min;
	      # @{ $local_proba_min[$pos_proba_min_0][0] } = @{ &clone($r_low_proba->[$pos_proba_min_0][$ilp]) };
	      $local_proba_min[$pos_proba_min_0][0] = $r_low_proba->[$pos_proba_min_0][$ilp];

	      $not_verbose or print "On lance EVAL line ".__LINE__.": $r_last_ind_seq_set->[$ilp]+1, $r_last_ind_seq_set->[$ilp+1], 0, $r_tempo_spacer_bounds_ind_l_this_pos, $pos_proba_min_0, $r_proba_ext, $r_border_words, $string_let_pos_p, 0, $regexp_motif_l, local_proba_min, $r_global_cpt_score_pos_proba_min, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice_l, local_ind_to_split\n";
	      # CAUTION: bool_not_recorded_interval set to 0 because we only want to test the motif: if bad evaluation, RECUP_LET_SORT ran some lines below, so we do not need
	      # to do it into EVAL (and would give wrong results as 2 RECUP_LE give a letter shift to high 
	      $bool_interess_res = &EVALUATION($r_last_ind_seq_set->[$ilp]+1, $r_last_ind_seq_set->[$ilp+1], 0, $r_tempo_spacer_bounds_ind_l_this_pos, $pos_proba_min_0, $r_proba_ext, $r_border_words, $string_let_pos_p, 0, \$regexp_motif_l, \@local_proba_min, $r_global_cpt_score_pos_proba_min, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice_l, \@local_ind_to_split);
	      
	    }
	    else{
	      $not_verbose or print "l'evaluation de $regexp_motif_l a deja du etre faite puisqu'aucune lettre n'a ete ajoutee!!!! line ".__LINE__."\n";
	     }

	    # we go on motif extension with interesting letter corresponding to the set only if we do not have already results for this set
	    if(not $bool_interess_res){
	      # **************************************************************************************************	    	      
	      $lettre = $rev_H_l[ $r_low_proba->[$pos_proba_min_0][$ilp][0] ]; # sorted tab, so letter with the lowest proba in first pos
	      # print "r_low_proba_pos_proba_min $r_low_proba->[$pos_proba_min_0][$ilp][0] doit correspondre a la lettre $lettre ln ".__LINE__."\n";	    
	      
	      # extension on the left
	      
	      # CAUTION nb_of_intern_l_added_to_boxes increased for futur calls &eval and &recup... DO NOT ADD other sort or treatments
	      if   ( $pos_proba_min_0 == 0 ){ $extended_border_words[0] = $lettre.$r_border_words->[ 0 ]; }
	      elsif( $pos_proba_min_0 == 2 ){ $extended_border_words[1] = $lettre.$r_border_words->[ 1 ]; $nb_of_intern_l_added_to_boxes++; }
	      # extension on the right
	      elsif( $pos_proba_min_0 == 1 ){ $extended_border_words[0] = $r_border_words->[ 0 ].$lettre; $nb_of_intern_l_added_to_boxes++; }
	      elsif( $pos_proba_min_0 == 3 ){ $extended_border_words[1] = $r_border_words->[ 1 ].$lettre; }
	      
	      
	      # + 1 corresponds to the letter which will be added thanks RECUP... useless as created regexp motif will be to long
	      if($nb_l_added_matrice + $initial_total_nb_boxes_letters + 1 >= $filtre_nb_l_added_matrice_maxi){
		
		($not_verbose)or print "EXTENSION TO LONG FOR futur regexp_motif, we do not need to get letters in RECUP , nb_l_added_matrice + initial_total_nb_boxes_letters ".($initial_total_nb_boxes_letters + $nb_l_added_matrice)." >= filtre_nb_l_added_matrice_maxi $filtre_nb_l_added_matrice_maxi\n";
		next;
	      }
	      
	      # $nb_l_added_matrice++;
	      my @extended_border_words_call = @extended_border_words;
	      
	      ($not_verbose)or print "On lance RECUP_LETTRE_SORT @bool sur les sequences indices ".($r_last_ind_seq_set->[$ilp]+1).", $r_last_ind_seq_set->[$ilp+1], dans sub SORT avec border_words @extended_border_words_call bool_not_recorded_interval_l ".(($#interv_int == -1)or($last_ind_seq_set_beg_l > $interv_int[$#interv_int]))." line ".__LINE__."\n";

	      # we go on with each letter (so seq set), we DO NOT TREAT fisrt best letter, then other
	      &RECUP_LETTRE_SORT($r_last_ind_seq_set->[$ilp]+1, $r_last_ind_seq_set->[$ilp+1], @bool , $bool_not_recorded_interval_l, \@extended_border_words_call, \@proba_ext, $bool_split_in_record_not_done, $string_let_pos_p, $nb_of_intern_l_added_to_boxes, $nb_l_added_matrice);

	    }
	  }
	  @extended_border_words = ();
	}
	
      }
    return $returned_val;
    if(! $not_verbose){
      
      print ' ' x length($string_let_pos_p);
      print "FIN SUB Evaluation @_\n";
    }
  } # end of sub EVALUATION

sub max_motif_length_and_P_MM0(\$\$$\$){
    my ($r_regexp_motif, $r_max_motif_length, $bool_want_MM0_P, $r_P_whole, $r_P_upstr) = @_;

    $$r_P_whole = $$r_P_upstr = 1;
    $$r_max_motif_length = 0;

    # normal letters
    while($$r_regexp_motif =~ /(?:^|[^\[])([ACGTacgt]+)(?:$|[^\]])/g){ 
	$$r_max_motif_length += length($1); 

	if($bool_want_MM0_P){
	    foreach(split //, $1){
		$$r_P_whole *= $MM0{$_}[$sp];
		$$r_P_upstr *= $MM0_upstr{$_}[$sp];
	    }
	}
    }
    while($$r_regexp_motif =~ /\[([ACGTacgt]+)\]/g){ 
	$$r_max_motif_length++; 
    
	if($bool_want_MM0_P){
	    my ($moy_P_whole, $moy_P_upstr) = (0)x2;
	    my $nb_P = 0;
	    foreach(split //, $1){
		$moy_P_whole += $MM0{$_}[$sp];
		$moy_P_upstr += $MM0_upstr{$_}[$sp];
		$nb_P++;
	    }
	    $$r_P_whole *= ($moy_P_whole/$nb_P);
	    $$r_P_upstr *= ($moy_P_upstr/$nb_P);
	}
    }
    # \w{1,12} or \w{2}
    while($$r_regexp_motif =~ /(\d+)\}/g)   { $$r_max_motif_length += $1; }
    # \w (one joker letter)
    while($$r_regexp_motif =~ /\\w(?:$|\w)/g){ $$r_max_motif_length++; }
    
    $not_verbose or print "max_motif_length $$r_max_motif_length, MM0 P_whole = $$r_P_whole, P_upstr $$r_P_upstr, line ".__LINE__."\n";
    # ?? VERIF DEBUG
}  # end of sub max_motif_length_and_P_MM0

# ****************************************************************************************************
# SUB FOR DEBUGGING **********************************************************************************

 sub verif_low_prob(\@$$){
   my ($r_low_proba, $pos_proba_min_0, $line) = @_;
   
   print "verif low_proba de pos_proba_min_0 $pos_proba_min_0 line $line\n";

   # liste de couples lettre, proba pour cette position
   for (my $j = 0; $j <= $#{$r_low_proba->[$pos_proba_min_0]}; $j++){
     
     print "pos_proba_min_0 $pos_proba_min_0, j $j, val @{$r_low_proba->[$pos_proba_min_0][$j]}\n";
   }
   if($#{$r_low_proba->[$pos_proba_min_0]} == -1){
     print "no record in r_low_proba $pos_proba_min_0 but ref defined\n";
   }

   print "end verif low_proba de pos_proba_min_0 $pos_proba_min_0\n";
 }

# ****************************************************************************************************
