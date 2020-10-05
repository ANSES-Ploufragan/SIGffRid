#!/usr/bin/perl
use warnings;
use strict;
use Env qw(SIGFFRID_DIR SIGFFRID_TEST_DIR); #  ANSES_MEMORY_LIMIT);
use Getopt::Long;
use Cwd 'abs_path';                             # to get absolute path from relative
# use Linux::MemInfo;                           # to get max RAM on machine
use Statistics::Basic qw(mean median stddev);
# use Statistics::Descriptive;
use POSIX qw(strftime);

use lib $SIGFFRID_DIR.'libperl';
use get_dir;    # get directory of a path
use short_name; # get name of file without path and extension
use run_capture;

# # **********************************************************************
$|++; # ensure display of print/errors in the order of appearing


# **********************************************************************
# For windows portability?
# see http://perldoc.perl.org/perlport.html for portability
# **********************************************************************
BEGIN {
    if ($^O =~ /^(MS)?Win/) {
        eval "use Win32::DriveInfo";
    }
}
# **********************************************************************

=head1  NAME

sigffrid_cmd_line.pl

=head1 DESCRIPTION

Run all scripts of SIGffRid program.

=head1  USAGE

=over

=item -gb_f1 <s>

Genbank file for first bacteria (must contain nucleotide sequence).

=item -gb_f2 <s>

Genbank file for second bacteria (must contain nucleotide sequence).

=item -id_bact1 <s>

ID of first bacteria (this id is in the file provided by MBGD database
file of ortholog comparisons.

=item -id_bact2 <s>

ID of second bacteria (this id is in the file provided by MBGD database
file of ortholog comparisons.

=item -MBGD_f <s>

MBGD file with ortholog relashionships between the genes of the two
bacteria.
Can be obtained by following SIGffRid instruction of README.md
file starting at this url
L<< MBGD_organism_selection|http://mbgd.genome.ad.jp/htbin/SelectOrganism.pl >>

=item [-filterout_16SrRNA_res]

To remove from results motifs related to ribosome binding site
(matching 16S rRNA genes and with a constant position compared to START pos).

=item [-cluster_mod <s>]

Default:no.

Other possible values:
'sge' to use SGE job scheduler or
'slurm' to use slurm.

=item [-res_dir <s>]

Default:directory of (first) file of reads directory where results are written.

=item [-nb_threads <i>]

Default:10 (caution, for 64 GB computer),  to use several threads when 
running assembler.

=item [-force]

To recompute an assembly already done (force erase previous results).

=item [-use_fct]

To create orhologs directories related to some functions (for specific sigma factors).

=item [-test]

To run only test for SIGffRid, do not take into account other arguments.

=item [-test_d_RBS_START]

To run only motif evaluation at the end to filter out RBS (test only on a toy
dataset, does not take into account other parameters).

=back

examples:


Will run for example following command lines:

perl sigffrid_cmd_line.pl -gb_f1 \
../ref_genomes/NC_002516.2_Pseudomonas_aeruginosa_PAO1_chr_complete_genome.gb
-gb_f2 ../ref_genomes/NC_009512.1_Pseudomonas_putida_F1_complete_genome.gb 
-id_bact1 pae -id_bact2 ppf -MBGD_f MBGD_file/cluster_PA_PPTU_RS.tab 
-res_dir sigffrid_res/


=cut

# **********************************************************************

my $verbose                        = 0;
my $b_force                        = 0;
my $b_test                         = 0;
my $b_test_d_RBS_START             = 0;
my $prog_tag                       = '[sigffrid_cmd_line.pl]';
my $upstream_seq_start_to_TSS      = 350;
my $overlap_upstr_seq_on_start     = 0;
my $bool_intergenic_for_evaluation = 0; # default 0
my $b_filterout_16SrRNA_res        = 0;
my $cluster_mod                    = 'no'; # can be 'sge' or 'slurm'
my $b_sge                          = 0;
my $b_slurm                        = 0;
my $display_width_before_column_ok = 85;

my $b_small_genome = 0;
my $nb_threads     = 10; # 64 GB computer!!!

my $start_date = undef;
my $end_date   = undef;

my $SIGffRID_scripts_dir = get_dir(abs_path(__FILE__));
# ok
# die "$prog_tag SIGffRID_scripts_dir:$SIGffRID_scripts_dir, line ".__LINE__."\n";

my @seeds = (
	'www',
	'wwww',
	'wnnw',
	'wnnnw',
	'wnw',
	'wwnnw',
	'wnnww',
	'wnwnw',
	'wwnw',
	'wnww'
);

# **********************************************************************
# variables
# **********************************************************************

my $gb_f1      = undef;
my $gb_f2      = undef;

my $id_bact1   = '';
my $id_bact2   = '';

my $MBGD_f     = undef;

my $sigffrid_res_dir = undef;

# previous 14-20 but some non sig70 promoters can have shorter spacers
# ideally 2 ranges: 14-20 and 7-13
my @min_spacer = (14,7); 
my @max_spacer = (20,13);
my $spacer_variation_1bact        = 1;
my $spacer_variation_between_bact = 1;

my $cmd    = undef;
my @res    = undef;
my $prefix = undef;

my $MM_order = 3; # >=3 <8

my @pids = ();

my $b_create_directories_for_fct = 0;

# ram slurm allocation MB
my $ram_Fextract_seq = 1000;
my $ram_rmes         = 2000;
my $ram_mbgd         = 500;
my $ram_Ftest_rech   = 4000;

# **********************************************************************


# **********************************************************************
# CHECK OPTIONS
# **********************************************************************
if((not defined $SIGFFRID_DIR)or($SIGFFRID_DIR =~ /^\s*$/))
{
    die "$prog_tag [Error] SIGFFRID_DIR      environment variable must be defined and not empty, line ".__LINE__."\n";
}
if((not defined $SIGFFRID_TEST_DIR)or($SIGFFRID_TEST_DIR =~ /^\s*$/))
{
    die "$prog_tag [Error] SIGFFRID_TEST_DIR environment variable must be defined and not empty, line ".__LINE__."\n";
}
if($SIGFFRID_DIR !~ /\/$/){      $SIGFFRID_DIR      .= '/' }
if($SIGFFRID_TEST_DIR !~ /\/$/){ $SIGFFRID_TEST_DIR .= '/' }

my $nbargmini = 1;
if(scalar(@ARGV) < $nbargmini){ 
  print "Bad number of arguments: ".scalar(@ARGV)." found, at least $nbargmini wanted\n";
  foreach(0..$#ARGV){
    print "$_: $ARGV[$_]\n";
  }
  die `perldoc $0`;
}
GetOptions(
    "gb_f1=s"                 => \$gb_f1,
    "gb_f2=s"                 => \$gb_f2,
    "id_bact1=s"              => \$id_bact1,
    "id_bact2=s"              => \$id_bact2,
    "MBGD_f=s"                => \$MBGD_f,
    "res_dir=s"               => \$sigffrid_res_dir,
    "nb_threads=i"            => \$nb_threads,
    "force"                   => sub { $b_force                      = 1  },
    "use_fct"                 => sub { $b_create_directories_for_fct = 1  },
    "filterout_16SrRNA_res"   => sub { $b_filterout_16SrRNA_res      = 1  },
    "cluster_mod=s"           => \$cluster_mod,
    "test"                    => sub { $b_test                       = 1  },
    "test_d_RBS_START"        => sub { $b_test_d_RBS_START           = 1  }
);
# **********************************************************************

foreach($cluster_mod)
{
    /sge/ and do{
	$b_sge = 1;
	$b_slurm = 0;	    	    
	last;
    };
    /slurm/ and do{
	$b_slurm = 1;
	$b_sge   = 0;
	last;
    };
    /no/ and do{
	$b_sge   = 0;
	$b_slurm = 0;	    
	last;
    };	
}

# *******************************************************
# subprog to deal with SGE processes
# *******************************************************
sub wait_processes_finished($)
{
    my ($r_pids) = @_;
    
    # wait for all process finished 
    while(scalar(@$r_pids) != 0)
    {
	my $b_test_f = 0;
	if($b_sge)     { $b_test_f = run_capture::is_run_capture_sge_finished(  $r_pids->[0], __LINE__) }
	elsif($b_slurm){ $b_test_f = run_capture::is_run_capture_slurm_finished($r_pids->[0], __LINE__) }
	else           { $b_test_f = run_capture::is_run_capture_finished(      $r_pids->[0], __LINE__) }
	
	if($b_test_f)
	{
	    shift @$r_pids;
	    next;
	}
	if(scalar(@$r_pids) > 0)
	{
	    print "$prog_tag waiting end of processes:".join(',', @$r_pids).", line ".__LINE__."\n";
	    sleep 1;
	}
    }
}

sub wait_processes_remaining_less_than($$)
{
    my($max_nb_processes_to_run, $r_pids) = @_;
    $max_nb_processes_to_run--;
    # wait for all process finished 
    while(scalar(@$r_pids) > $max_nb_processes_to_run)
    {
	my $b_test_f = 0;
	if($b_sge)     { $b_test_f = run_capture::is_run_capture_sge_finished(  $r_pids->[0], __LINE__) }
	elsif($b_slurm){ $b_test_f = run_capture::is_run_capture_slurm_finished($r_pids->[0], __LINE__) }
	else           { $b_test_f = run_capture::is_run_capture_finished(      $r_pids->[0], __LINE__) }
	
	if($b_test_f)
	{
	    shift @$r_pids;
	    next;
	}
	scalar(@$r_pids) and sleep 1;
    }
}
# *******************************************************


# *******************************************************
# getCosmidSequence_gb [f] returns the complete cosmid nucleotid sequence from EMBL format file [f]. 
# Resulting sequence is one big lowercase nucleotid string.
# *******************************************************
sub getCosmidSequence_gb($)
{
  my $fileName = shift @_;
  my $result = ""; #Result buffer.
  my ($currentLine,$sequenceBlock);
  
  open(GBFILE, '<', $fileName) or die "$prog_tag [Error] Can't open GB GBK file \"$fileName\": $!, line ".__LINE__."\n"; #GB GBK format file.
  if($fileName =~ /\.(?:gbk|gb)$/){
      while($currentLine = <GBFILE>)
      {
	  if(($sequenceBlock) = ($currentLine =~ /^\s+\d+([\satgc]+)$/)) #Current line is a DNA sequence line.
	  {
	      $sequenceBlock =~ s/\s+//g; #Eat up whitespace.
	      $result .= $sequenceBlock ;#Add the DNA block to final sequence.	      
	  }
      }
      
  }
  else{
      die "$prog_tag [Error] You must give a filename with extensions .gb or .gbk !\n";
  }
  close GBFILE;
  return $result;
}
# *******************************************************

# *******************************************************
# subprog to filter out motif with 16S RNA sequence
# if no filter is required, file name does not change, otherwise
# a new file with name adding _RBSoufiltered before .txt extension 
# is created
# NEED TO ADD ANALYSE OF DISTANCE TO START (CONSTANT FOR RBS)
# *******************************************************
sub filterout_16SrRNA_res($$$$)
{
    use SeqUtilities qw(complement_seq);
    my($grep_sortie_f, 
       $gb_f,
       $user_bact,
       $dir_of_fasta_for_motifs
	) = @_;
	my $n_after_product = 2;
    my $start_16S = undef;
    my $end_16S   = undef;	    
    my $grep_sortie_filtered_f = $grep_sortie_f;
    $grep_sortie_filtered_f =~ s/.txt$/_RBSfilteredout.txt/;
    if($grep_sortie_f eq $grep_sortie_filtered_f)
    {
	die "$prog_tag [Error] Cannot deduce grep_sortie_filtered_f from $grep_sortie_f, probably does not have .txt extension, line ".__LINE__."\n";
    }

    my $max_stderr_for_RBS = 6;
    
    # get fasta sequence from whole genome/gb file
    my $seq_16S = getCosmidSequence_gb($gb_f);
    
    # ----------------------------------------------------		    
    # search first coordinates of 16S sRNA in GenBank file to extract nt sequence
    $cmd = "egrep -i -A $n_after_product 'product=\"1[68]S ribosomal RNA|product=\"Ribosomal RNA 1[68]S' $gb_f";
    my @res = ();
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    $prefix = 'egrep_16_18S_in_gb'; 
    if($b_sge)     {  run_capture::run_capture_get_res_sge_prefix_wait_finished(  $cmd, __LINE__, $prefix, \@res)     ; }
    elsif($b_slurm){  run_capture::run_capture_get_res_slurm_prefix_wait_finished($cmd, __LINE__, $prefix, 200, \@res); }
    else           {  run_capture::run_capture_get_res(                           $cmd, __LINE__, \@res)              ; } 
    print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Search for 16/18S sRNA in GenBank file, ran").": ok\n");

	for my $i(1..$#res)
	{
	    if($res[$i] =~ /^[\s\t]+gene[\s\t]+(\d+)\.\.(\d+)\s*$/)
	    {
		$start_16S = $1;
		$end_16S   = $2;
		$seq_16S   = substr($seq_16S, $start_16S - 1, $end_16S - $start_16S + 1);
		last;
	    }
	    elsif($res[$i] =~ /^[\s\t]+gene[\s\t]+complement\((\d+)\.\.(\d+)\)\s*$/)
	    {
		$start_16S = $1;
		$end_16S   = $2;
		$seq_16S   = complement_seq(substr($seq_16S, $start_16S - 1, $end_16S - $start_16S + 1));
		last;
	    }
	}

	if(not defined $start_16S)
	{
	    warn "$prog_tag [filterout_16SrRNA_res][WARN] do not find product=\"1[68]S ribosomal RNA\" in $gb_f file, maybe modifiy grep to get not only first following line but more to get line with coordinates, line ".__LINE__."\nMaybe increase the number of lines analysed following 'product' to get 'gene' line with coordinates (currently $n_after_product lines), line ".__LINE__."\n";
	}
	else
	{
	    print "$prog_tag [filterout_16SrRNA_res] 16/18S sequence obtained from $gb_f:$seq_16S\n";
	}
    # ----------------------------------------------------
    
    # ----------------------------------------------------
    # compare 16S seq with motifs found in grep_sortie_f file
    # retain only lines with motifs not related to RBS (means these are not shine delgarno boxes)
    # he RBS typically lies about 7 nucleotides upstream from the start codon, and is AGGAGG
    # (data from http://parts.igem.org/Help:Ribosome_Binding_Site)
    open(ORIGREPOUT,'<',$grep_sortie_f)or die "$prog_tag [Error] Cannot open $grep_sortie_f file:$!, line ".__LINE__."\n";
    open(NEWGREPOUT,'>',$grep_sortie_filtered_f)or die "$prog_tag [Error] Cannot create $grep_sortie_filtered_f file:$!, line ".__LINE__."\n";    
    while(my $grepline = <ORIGREPOUT>)
    {
	# lines results from "grep 'MOTIF' ${sigffrid_res_dir}SIGffRid_orthologs/*/align*${id_bact}* | sort -nrk 4" that gives
	# SIGffRid_orthologs/align_14_20_1_1_www_wnnnw_gm04609.txt:MOTIF aata\w{12,15}agggg, R: 0.50  (>= 0.49  ), LRT is 11.99  (>= 3.84 , a = 0.05 ) nb seq 9, 30 in promot
	# motif is in second column if we separate by ' '
	my ($pathfile,$motif,@trash) = split /\s+/, $grepline;
	$motif =~ s/,$//;
	my $motif_right_part         = '';
	my $motif_right_part_revcomp = '';	
	if($motif =~ /([acgt]{4,})$/)
	{
	    $motif_right_part = $1;
	    # 16S can be ontained reverse complemented, not to reverse complement it,
	    # we reverse complement th motif to search (shorter)
	    $motif_right_part_revcomp = complement_seq($motif_right_part);
	}
	
	if(
	    ($motif_right_part ne '')and
	    # $b_test_d_RBS_START or
	   ($seq_16S =~ /(?:$motif_right_part|$motif_right_part_revcomp)/))
	{

	    # --------------------------------------------
	    # check distribution of positions of the motif

	    # get headers with motif position
	    my $motif_f = $dir_of_fasta_for_motifs.'motif'.$user_bact.'_'.$motif.'.txt';
	    # we get only part of header with the position relative to traduction start (last value on the line)
	    $cmd = "grep '>' '$motif_f'";
	    my @grep_res   = ();	    
	    my @pos_distri = ();
	    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
	    if($b_sge)     {  run_capture::run_capture_get_res_sge_prefix_wait_finished(  $cmd, __LINE__, $prefix, \@grep_res);      }
	    elsif($b_slurm){  run_capture::run_capture_get_res_slurm_prefix_wait_finished($cmd, __LINE__, $prefix, 200, \@grep_res); }
	    else           {  run_capture::run_capture_get_res(                           $cmd, __LINE__, \@grep_res);               }
	    foreach(@grep_res)
	    {
		/([\d\-]+)$/ and push @pos_distri, $1;
	    }
	    chomp(@pos_distri);
	    # numerical sort
	    print "$prog_tag pos_distri:".join(',',@pos_distri).", line ".__LINE__."\n";	    
	    @pos_distri = sort { $a <=> $b } @pos_distri;
	    print "$prog_tag pos_distri:".join(',',@pos_distri).", line ".__LINE__."\n";
	    # basic stats on pos
	    my $mediane    = median(\@pos_distri);
	    my $mean       = mean(\@pos_distri);
	    # for quantiles
	    my $Xn = $#pos_distri+1;
	    my $q0 = int(.1 * ($Xn + 1)) - 1;	    
	    my $q1 = int(.25 * ($Xn + 1)) - 1;
	    my $q2 = int(.5 * ($Xn + 1)) - 1;
	    my $q3 = int(.75 * ($Xn + 1)) - 1;
	    my $q4 = int(.9 * ($Xn + 1)) - 1;
	    print "$prog_tag q0,q1,q2,q3,q4:".join(',',$q0,$q1,$q2,$q3,$q4).", line ".__LINE__."\n";

	    my @pos_distri_50pc_lowest_vals = @pos_distri[$q2..$#pos_distri];
	    my $squareroot = stddev(\@pos_distri_50pc_lowest_vals);; # stddev(\@pos_distri);
	    
	    print join("\t", "$prog_tag [filterout_16SrRNA_res] $motif:", 
		       'median:'.$mediane,
		       'stderr 50% of val near START:'.$squareroot,
		       'mean:'.$mean,
		       # "quartiles $x10,$x25,$x50,$x75,$x90 (10,25,50,75,90)\n"
		       "freq_distri (10,25,50,75,90):".join(',',@pos_distri[$q0,$q1,$q2,$q3,$q4])."\n"		       
		);
	    
	    print "$prog_tag ran\n";
	    # --------------------------------------------
	    if($squareroot <= $max_stderr_for_RBS)
	    {
		# if too constant?, we remove because probable RBS, otherwise we keep
		print "$prog_tag [filterout_16SrRNA_res] remove $motif because found in ribosomal RNA\n";
	    }
	    else
	    {
		print "$motif included in 16/18S sequence portion but std dev is too large to have a RBS (Shine Delgarno motif), line ".__LINE__."\n";
		print NEWGREPOUT $grepline;
	    }
	}
	else
	{
	    print "$motif not included in 16/18S sequence, line ".__LINE__."\n";
	    print NEWGREPOUT $grepline;
	}
    }
    close(ORIGREPOUT);
    close(NEWGREPOUT);    
    # ----------------------------------------------------
	print "$prog_tag [filterout_16SrRNA_res] $grep_sortie_filtered_f file created\n";
    return $grep_sortie_filtered_f;
}
# ********************************************************

# --------------------------------------------------------
# program test
# --------------------------------------------------------
# my $sigffrid_test_dir = $SIGFFRID_DIR.'test_SIGffRid/';
# my $sigffrid_test_dir = '/media/sde_res2/test_SIGffRid/';
my $sigffrid_test_dir = $SIGFFRID_TEST_DIR.'test_SIGffRid/';

if($b_test_d_RBS_START)
{
    my $RBS_START_test_dir = $sigffrid_test_dir.'test_d_RBS_START/';
    
    my @motifs = ( # for sah Staphylococcus aureus
#		   "motifsah_tgac\w{18,20}tataat.txt", # sigma file
		   'tgac\w{18,20}tataat', # sigma motif sah
		   'tgttataat\w{12,14}a\wt', # sigma sah
		   't\wta\w{17,18}gttataat', # sigma sah
		   't\wt\w{17,19}aaggag', # RBS ? sah
		   'g\wa\w{18,20}aaggag', # RBS ? sah
		   't\wg\w{13,15}aaggag', # RBS ? sah
		  # for mtu (Mycobacterium tuberculosis)
		  'cgacgat\w{17,19}agga', # sigma mtu?
		  'gtta\w{8,8}taac', # other TFBS mtu
		   'acgatg\w{18,20}agga',  # sigma or RBS mtu ?
		   'a\wc\wa\w{5,10}aggagc' # match to RBS seq mtu
	);

	
    # for mtu
    my $user_bact = 'mtu';
    my $grep_sortie_mtu = $RBS_START_test_dir.'grep_sortie_mtu.txt';
	# '/media/sde_res2/test_SIGffRid/sigffrid_test_res/SIGffRid_orthologs/align_14_20_1_1_www_wwww_mtu.txt:9561:MOTIF cgacgat\w{17,19}agga, R: 0.62  (>= 0.48  ), LRT is 28.43  (>= 3.84 , a = 0.05 ) nb seq 8, 31 in promot',
	# '/media/sde_res2/test_SIGffRid/sigffrid_test_res/SIGffRid_orthologs/align_7_13_1_1_www_wwww_mtu.txt:23539:MOTIF gtta\w{8,8}taac, R: 0.56  (>= 0.48  ), LRT is 13.31  (>= 3.84 , a = 0.05 ) nb seq 9, 19 in promot',
    # '/media/sde_res2/test_SIGffRid/sigffrid_test_res/SIGffRid_orthologs/align_14_20_1_1_www_wwnw_mtu.txt:22375:MOTIF acgatg\w{18,20}agga, R: 0.65  (>= 0.48  ), LRT is 55.01  (>= 3.84 , a = 0.05 ) nb seq 8, 60 in promot'
    -e $grep_sortie_mtu or die "$prog_tag [test_d_RBS_START][Error] $grep_sortie_mtu file does not exist, please check test data, line ".__LINE__."\n";        
    my $gb_f = $sigffrid_test_dir.'ref_genomes/AL123456.3_Mycobacterium_tuberculosis_H37Rv_complete_genome.gb';
    -e $gb_f or die "$prog_tag [test_d_RBS_START][Error] $gb_f file does not exist, please check test data, line ".__LINE__."\n";
    my $dir_of_fasta_for_motifs = $RBS_START_test_dir;
    print "$prog_tag [test_d_RBS_START] run filterout_16SrRNA_res for $user_bact\n";    
    $grep_sortie_mtu = filterout_16SrRNA_res(
	$grep_sortie_mtu,
	$gb_f,
	$user_bact,
	$dir_of_fasta_for_motifs
	);
    
    # for sah
    $user_bact = 'sah';    
    my $grep_sortie_sah = $RBS_START_test_dir.'grep_sortie_sah.txt';
    -e $grep_sortie_sah or die "$prog_tag [test_d_RBS_START][Error] $grep_sortie_sah file does not exist, please check test data, line ".__LINE__."\n";    
    $gb_f = '/media/sde_res2/test_SIGffRid/test_Staphylococcus_argenteus_vs_Staphylococcus_capitis/Staphylococcus_aureus_subsp_aureus_JH1_NC_009632.1_complete_genome.gb';
    -e $gb_f or die "$prog_tag [test_d_RBS_START][Error] $gb_f file does not exist, please check test data, line ".__LINE__."\n";    
    $dir_of_fasta_for_motifs = $RBS_START_test_dir;
    print "$prog_tag [test_d_RBS_START] run filterout_16SrRNA_res for $user_bact\n";
    $grep_sortie_sah = filterout_16SrRNA_res(
	$grep_sortie_sah,
	$gb_f,
	$user_bact,	
	$dir_of_fasta_for_motifs
	);
    
    exit;
}
if($b_test)
{
    $sigffrid_res_dir = $sigffrid_test_dir.'sigffrid_test_res/';
    
    $id_bact1 = 'mtu';
    $id_bact2 = 'mmi';
    $gb_f1    = $sigffrid_test_dir.'ref_genomes/AL123456.3_Mycobacterium_tuberculosis_H37Rv_complete_genome.gb';
    $gb_f2    = $sigffrid_test_dir.'ref_genomes/NC_010612.1_Mycobacterium_marinum_M_complete_genome.gb';
    $MBGD_f   = $sigffrid_test_dir.'MBGD_file/cluster_mtu_mmi.tab';
    
    foreach(	$sigffrid_test_dir,
		$sigffrid_res_dir,
		$gb_f1,
		$gb_f2,
		$MBGD_f)
    {
        -e $_ or die "$prog_tag [Error] $_ file/dir does not exist for TEST, have you moved it?, line ".__LINE__."\n";
    }
    
    my $short_name_gb_f1 = short_name($gb_f1);
    my $short_name_gb_f2 = short_name($gb_f2);
    
    my $prefix = 'test_sigffrid_cmd_line_';
    my $cmd = join(' ', "perl ${SIGFFRID_DIR}sigffrid_cmd_line.pl",
		   "-gb_f1 $gb_f1",
		   "-gb_f2 $gb_f2",
		   "-id_bact1 $id_bact1",
		   "-id_bact2 $id_bact2",
		   "-MBGD_f $MBGD_f",
		   "-res_dir $sigffrid_res_dir",
		   "-nb_threads $nb_threads");
    defined $cluster_mod             and $cmd .= " -cluster_mod $cluster_mod ";
    defined $b_filterout_16SrRNA_res and $cmd .= " -filterout_16SrRNA_res ";
    
    print "$prog_tag [TEST] cmd:$cmd, line ".__LINE__."\n";
    # if($b_sge)     { run_capture::run_capture_sge_prefix_wait_finished(  $cmd, __LINE__, $prefix); }
    # elsif($b_slurm){ run_capture::run_capture_slurm_prefix_wait_finished($cmd, __LINE__, $prefix); }
    # else           { run_capture::run_capture_prefix_wait_finished(      $cmd, __LINE__, $prefix); }
    `$cmd`;
    print "$prog_tag [TEST] ran\n";
    exit;
}
# --------------------------------------------------------
else
{
    # get start date, something like "Thu Oct 13 04:54:34 1994"
    $start_date = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print "$prog_tag start:$start_date\n";

    # normal run
    (defined $gb_f1)            or die "$prog_tag [Error] gb_f1 file/directory not defined, line ".__LINE__."\n";
    (defined $gb_f2)            or die "$prog_tag [Error] gb_f2 file/directory not defined, line ".__LINE__."\n";
    (defined $id_bact1)         or die "$prog_tag [Error] id_bact1 file/directory not defined, line ".__LINE__."\n";
    (defined $id_bact2)         or die "$prog_tag [Error] id_bact2 file/directory not defined, line ".__LINE__."\n";
    (defined $MBGD_f)           or die "$prog_tag [Error] MBGD_f file/directory not defined, line ".__LINE__."\n";
    (defined $sigffrid_res_dir) or die "$prog_tag [Error] res_dir file/directory not defined, line ".__LINE__."\n";
    
    # check if defined files/directories exist
    foreach(grep { defined $_ } 	$gb_f1,
	    $gb_f2,
	    $MBGD_f,
	    $sigffrid_res_dir)
    {
        -e $_ or die "$prog_tag [Error] $_ file/directory does not exist, line ".__LINE__."\n";
    }
    
    my $short_name_gb_f1 = short_name($gb_f1);
    my $short_name_gb_f2 = short_name($gb_f2);
    
    my $orthologs_dir = $sigffrid_res_dir.'SIGffRid_orthologs/';
    my $final_res_dir = $sigffrid_res_dir.'SIGffRid_results/';
    
    # --------------------------------------------------------
    # formatting sequences for bact1 and bact2 in parallel
    # bact1
    $cmd = join(' ', 	"perl ${SIGffRID_scripts_dir}FextractSeqFeaturesFromGbUsingBioPerl.pl",
                '-id_bact', $id_bact1,
                '-gb_f'   , $gb_f1,
                '-res_dir', $sigffrid_res_dir);
    $prefix = 'FextractSeqFeaturesFromGbUsingBioPerl_'.$id_bact1.'_'.$short_name_gb_f1;
    
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    if($b_sge)     {  push @pids, run_capture::run_capture_sge_prefix(  $cmd, __LINE__, $prefix);                    }
    elsif($b_slurm){  push @pids, run_capture::run_capture_slurm_prefix($cmd, __LINE__, $prefix, $ram_Fextract_seq); }
    else           {  push @pids, run_capture::run_capture_prefix(      $cmd, __LINE__, $prefix);                    }    
    
    # bact2
    $cmd = join(' ', 	"perl ${SIGffRID_scripts_dir}FextractSeqFeaturesFromGbUsingBioPerl.pl",
                '-id_bact', $id_bact2,
                '-gb_f'   , $gb_f2,
                '-res_dir', $sigffrid_res_dir);
    $prefix = 'FextractSeqFeaturesFromGbUsingBioPerl_'.$id_bact2.'_'.$short_name_gb_f2;
    
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    if($b_sge)     {  push @pids, run_capture::run_capture_sge_prefix(  $cmd, __LINE__, $prefix);                    }
    elsif($b_slurm){  push @pids, run_capture::run_capture_slurm_prefix($cmd, __LINE__, $prefix, $ram_Fextract_seq); }
    else           {  push @pids, run_capture::run_capture_prefix(      $cmd, __LINE__, $prefix);                    } 
    
    wait_processes_finished(\@pids);
    print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Sequences extraction for $id_bact1 and $id_bact2 Search for 16/18S sRNA in GenBank file, ran").": ok\n");
    
    # get name of merged upstream seq files normally created by previous
    # script, used by rmes
    #my $merged_upstr_f1 = $sigffrid_res_dir.'not_fasta_-350_0_annot_overlapped_seq_sens01_'.$id_bact1.'.txt';
    #my $merged_upstr_f2 = $sigffrid_res_dir.'not_fasta_-350_0_annot_overlapped_seq_sens01_'.$id_bact2.'.txt';
    my $merged_upstr_f1 = $sigffrid_res_dir.'not_fasta_-'.$upstream_seq_start_to_TSS.'_'.$overlap_upstr_seq_on_start.'_overlapped_seq_sens01_'.$id_bact1.'.txt';
    my $merged_upstr_f2 = $sigffrid_res_dir.'not_fasta_-'.$upstream_seq_start_to_TSS.'_'.$overlap_upstr_seq_on_start.'_overlapped_seq_sens01_'.$id_bact2.'.txt';
    
    # get name of upstream seq files normally created by previous script
    # promoters of intergenic regions
    my $upstr_f1 = $sigffrid_res_dir.'not_fasta_promot_'.$id_bact1.'.txt';
    my $upstr_f2 = $sigffrid_res_dir.'not_fasta_promot_'.$id_bact2.'.txt';
    
    # verif that files exist
    foreach($merged_upstr_f1,
	    $merged_upstr_f2,
	    $upstr_f1,
	    $upstr_f2)
    {
	-e $_ or die "$prog_tag [Error] $_ file does not exist, maybe $prefix failed, line ".__LINE__."\n";
    }
    
    
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    # markov process for bact1 and bact2 in parallel

    # expected out files
    my $markov_expected_out_f1 = undef;
    my $markov_expected_out_f2 = undef;
    if($gb_f1 =~ /([^\/]+)\.(?:embl|gbk|gb|txt|fasta|fa|GBK|fsa|fa\.gz)$/)
    {
        $markov_expected_out_f1 = join('', $sigffrid_res_dir,'markov_mod_order',$MM_order,'_',$1,'.txt');	
    }
    else
    {
        die "$prog_tag [Error] Unrecognized gb file: $gb_f1, line ".__LINE__."\n";
    }
    if($gb_f2 =~ /([^\/]+)\.(?:embl|gbk|gb|txt|fasta|fa|GBK|fsa|fa\.gz)$/)
    {
        $markov_expected_out_f2 = join('', $sigffrid_res_dir,'markov_mod_order',$MM_order,'_',$1,'.txt');	
    }
    else
    {
        die "$prog_tag [Error] Unrecognized gb file: $gb_f2, line ".__LINE__."\n";
    }
    
    # bact1
    $cmd = join(' ', 	"perl ${SIGffRID_scripts_dir}FmarkovBackgroundSeveralBacteria4.pl",
                $MM_order,
                $gb_f1,
                $sigffrid_res_dir,
                0);
    $prefix = 'FmarkovBackgroundSeveralBacteria4_'.$id_bact1;
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    if($b_sge)     {  push @pids, run_capture::run_capture_sge_prefix(  $cmd, __LINE__, $prefix);       }
    elsif($b_slurm){  push @pids, run_capture::run_capture_slurm_prefix($cmd, __LINE__, $prefix, $ram_rmes); }
    else           {  push @pids, run_capture::run_capture_prefix(      $cmd, __LINE__, $prefix);       } 
    
    
    # bact2
    $cmd = join(' ', 	"perl ${SIGffRID_scripts_dir}FmarkovBackgroundSeveralBacteria4.pl",
                $MM_order,
                $gb_f2,
                $sigffrid_res_dir,
                0);
    $prefix = 'FmarkovBackgroundSeveralBacteria4_'.$id_bact2;
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    if($b_sge)     {  push @pids, run_capture::run_capture_sge_prefix(  $cmd, __LINE__, $prefix);            }
    elsif($b_slurm){  push @pids, run_capture::run_capture_slurm_prefix($cmd, __LINE__, $prefix, $ram_rmes); }
    else           {  push @pids, run_capture::run_capture_prefix(      $cmd, __LINE__, $prefix);            } 
    
    wait_processes_finished(\@pids);
    print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Markov models computation for $id_bact1 and $id_bact2, ran").": ok\n");
    
    # verif that output files exist, check syntax
    foreach($markov_expected_out_f1, $markov_expected_out_f1)
    {
    	-e $_ or die "$prog_tag [Error] $_ file does not exist, maybe $prefix failed, line ".__LINE__."\n";
    }
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    # rmes process for bact1 and bact2 in parallel
    my $RMES_res_dir = $sigffrid_res_dir.'RMES/';
    -e $RMES_res_dir or mkdir $RMES_res_dir;

    # expected RMES output files
    my $RMES_out_expected_f1 = "${RMES_res_dir}mots_plus_comp_sans_palind_${id_bact1}_rmes_8_3.txt";
    my $RMES_out_expected_f2 = "${RMES_res_dir}mots_plus_comp_sans_palind_${id_bact2}_rmes_8_3.txt";

    # bact1
    $cmd = join(' ', 	"perl ${SIGffRID_scripts_dir}RMES/Frmes_fam_select_gauss_poiss.pl",
                $id_bact1, # suffix for created r'mes files
                $merged_upstr_f1, # name of the sequence file
                # TO DO: CHANGE LOCATION OF RMES RES
                $RMES_res_dir);
    $prefix = 'Frmes_fam_select_gauss_poiss_'.$id_bact1;
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    if($b_sge)     {  push @pids, run_capture::run_capture_sge_prefix(  $cmd, __LINE__, $prefix);             }
    elsif($b_slurm){  push @pids, run_capture::run_capture_slurm_prefix($cmd, __LINE__, $prefix,  $ram_rmes); }
    else           {  push @pids, run_capture::run_capture_prefix(      $cmd, __LINE__, $prefix);             } 

    # bact2
    $cmd = join(' ', 	"perl ${SIGffRID_scripts_dir}RMES/Frmes_fam_select_gauss_poiss.pl",
                $id_bact2, # suffix for created r'mes files
                $merged_upstr_f2, # name of the sequence file
                # TO DO: CHANGE LOCATION OF RMES RES
                $RMES_res_dir);
    $prefix = 'Frmes_fam_select_gauss_poiss_'.$id_bact2;
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    if($b_sge)     {  push @pids, run_capture::run_capture_sge_prefix(  $cmd, __LINE__, $prefix);            }
    elsif($b_slurm){  push @pids, run_capture::run_capture_slurm_prefix($cmd, __LINE__, $prefix, $ram_rmes); }
    else           {  push @pids, run_capture::run_capture_prefix(      $cmd, __LINE__, $prefix);            } 
    
    wait_processes_finished(\@pids);
    print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag R'MES computation for $id_bact1 and $id_bact2, ran").": ok\n");

    print("RMES_out_expected_f1:$RMES_out_expected_f1, line ".__LINE__."\n");

    # verif that files exist, check syntax
    foreach($RMES_out_expected_f1, $RMES_out_expected_f2)
    {
        -e $_ or die("$prog_tag [Error] $_ file does not exist, maybe $prefix failed, line ".__LINE__."\n");
    }
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    # data extraction from orthologous files
    $cmd = join(' ',	"perl ${SIGffRID_scripts_dir}Ftraitement_fichier_sortie_MBGD.pl",
		$MBGD_f,
		$sigffrid_res_dir,
		$id_bact1,
		$id_bact2);
    $prefix = 'Ftraitement_fichier_sortie_MBGD_'.short_name($MBGD_f);
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    if($b_sge)     {  run_capture::run_capture_sge_prefix_wait_finished(  $cmd, __LINE__, $prefix);            }
    elsif($b_slurm){  run_capture::run_capture_slurm_prefix_wait_finished($cmd, __LINE__, $prefix, $ram_mbgd); }
    else           {  run_capture::run_capture_prefix_wait_finished(      $cmd, __LINE__, $prefix);            } 
    print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Orthologous deduction for $id_bact1 and $id_bact2, ran").": ok\n");            
    
    # file created by previous script, we check it exists
    my $file_with_gene_ids = $sigffrid_res_dir.join('_', 	
						    'ortho',
						    $id_bact1,
						    $id_bact2,
						    'fct_match.txt');
    
    -e $file_with_gene_ids or die "$prog_tag [Error] $file_with_gene_ids file does not exist, Ftraitement_fichier_sortie_MBGD.pl failed? , please check, line ".__LINE__."\n";
    
    -e $orthologs_dir or mkdir $orthologs_dir;
    
    $cmd = join(' ',	"perl ${SIGffRID_scripts_dir}Fprend_seq_prom_ortho_pairs_fct.pl",
		$file_with_gene_ids, # <file with gene ID (blast_donne_orth... or ortho...)>
		$id_bact1,
		$id_bact2,
		$upstr_f1, # <file containing upstream sequences of the first given organism>
		$upstr_f2, # <file containing upstream sequences of the second given organism>
		$b_create_directories_for_fct,
		$sigffrid_res_dir);
    $prefix = 'Fprend_seq_prom_ortho_pairs_fct'.short_name($file_with_gene_ids);
    print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
    if($b_sge)     {  run_capture::run_capture_sge_prefix_wait_finished(  $cmd, __LINE__, $prefix);      }
    elsif($b_slurm){  run_capture::run_capture_slurm_prefix_wait_finished($cmd, __LINE__, $prefix, $ram_mbgd); }
    else           {  run_capture::run_capture_prefix_wait_finished(      $cmd, __LINE__, $prefix);      } 
    print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Extraction of orthologous promoters for $id_bact1 and $id_bact2, ran").": ok\n");            

    # # verif that files exist, check syntax, not possible, depends on treatment, gene names...
    # foreach()
    # {
    # 	-e $_ or die "$prog_tag [Error] $_ file does not exist, maybe $prefix failed, line ".__LINE__."\n";
    # }
    
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    # real sigffrid analysis
    my $rmes_res_dir = $sigffrid_res_dir.'RMES/';
    -e $rmes_res_dir or mkdir $rmes_res_dir;
    my $second_prefix_of_promoter_f = "-".$upstream_seq_start_to_TSS.'_'.$overlap_upstr_seq_on_start.'_';
    print "$prog_tag second_prefix_of_promoter_f:$second_prefix_of_promoter_f, line ".__LINE__."\n";

    foreach my $i_spacer_range(0..$#min_spacer)
    {
	my $min_spacer = $min_spacer[$i_spacer_range];
	my $max_spacer = $max_spacer[$i_spacer_range];

	# this boolean is inititialized each time spacer parameters change
	my $bool_www_treated = 0;
	foreach my $seed1_i(0..$#seeds)
	{
	    my $seed1 = $seeds[$seed1_i];
	    foreach my $seed2_i($seed1_i +1..$#seeds)
	    {
		my $seed2 = $seeds[$seed2_i];

		# if pair of plain seeds already treated, no need to run it again
		if($bool_www_treated and (($seed1 =~ /w{3}/)and($seed2 =~ /w{3}/)))
		{
		    next;
		}
		$cmd = join(' ', "perl ${SIGffRID_scripts_dir}Ftest_rech_motif_fic_diff_ds_seq_promot_averees_vari_spacer_triTwoBoxes64multi.pl",
			    $id_bact2,
			    $id_bact1,
			    $orthologs_dir,
			    "${RMES_res_dir}mots_plus_comp_sans_palind_${id_bact1}_rmes_8_3.txt",
			    "${sigffrid_res_dir}markov_mod_order${MM_order}_${short_name_gb_f1}.txt",
			    "${sigffrid_res_dir}markov_mod_order${MM_order}_${short_name_gb_f2}.txt",
			    $min_spacer,
			    $max_spacer,
			    $spacer_variation_1bact,
			    $spacer_variation_between_bact,
			    $bool_intergenic_for_evaluation, # 0,
			    $sigffrid_res_dir,
			    $seed1, # www
			    $seed2, # wnww
			    0, # $bool_www_treated, # boolean to tell if www www combination was already treated
			    'not_fasta_',                # first prefix of promoter files
			    $second_prefix_of_promoter_f # second prefix of promoter files
		    );
		$prefix = join('_',	'Ftest_rech_motif',
			       $id_bact2,
			       $id_bact1,
			       $seed1, # www
			       $seed2);
		
		wait_processes_remaining_less_than($nb_threads, \@pids);
		print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
		if($b_sge)     { push @pids, run_capture::run_capture_sge_prefix(  $cmd, __LINE__, $prefix);                  }
		elsif($b_slurm){ push @pids, run_capture::run_capture_slurm_prefix($cmd, __LINE__, $prefix, $ram_Ftest_rech); }
		else           { push @pids, run_capture::run_capture_prefix(      $cmd, __LINE__, $prefix);                  } 
		print "$prog_tag process $pids[$#pids] added, line ".__LINE__."\n";
		
		if(($seed1 =~ /\w{3}/)and($seed2 =~ /\w{3}/))
		{
		    $bool_www_treated = 1;
		}
	    }
	}
    }
    print "$prog_tag wait end of processes:\n".join(',',@pids).", line ".__LINE__."\n";    
    wait_processes_finished(\@pids);
    print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Searchs of conserved motifs and their extension for $id_bact1 and $id_bact2, ran").": ok\n");     
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    # grep results for final analysis files
    foreach my $id_bact($id_bact1, $id_bact2)
    { 
	
	# this part we use to make the align_id.txt files that contains the final motifs in a file
	my $final_results_dir = $sigffrid_res_dir.'SIGffRid_results/FINAL_RESULTS/';
	-e $final_results_dir or mkdir($final_results_dir);
	my $grep_sortie= $final_results_dir.'result_grep_'.$id_bact.'.txt';
	(-e $grep_sortie) and unlink($grep_sortie);

	if($b_create_directories_for_fct)
	{
	    $cmd = join(' ', 	"grep 'MOTIF' ${sigffrid_res_dir}SIGffRid_orthologs/*/align*${id_bact}* |",
			# "tr '.' ',' | ", # only for FR
			'sort -nrk 4', # sort by R
			"> $grep_sortie");
	}
	else
	{
	    $cmd = join(' ', 	"grep -n 'MOTIF' ${sigffrid_res_dir}SIGffRid_orthologs/align*${id_bact}* |",
			# "tr '.' ',' | ", # only for FR
			'sort -nrk 4', # sort by R
			"> $grep_sortie");
	}
	$prefix = 'grep_MOTIF_'.short_name($file_with_gene_ids);
	print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
	if($b_sge)     {  run_capture::run_capture_sge_prefix_wait_finished($cmd, __LINE__, $prefix);              }
	elsif($b_slurm){  run_capture::run_capture_slurm_prefix_wait_finished($cmd, __LINE__, $prefix, $ram_mbgd); }
	else           {  run_capture::run_capture_prefix_wait_finished($cmd, __LINE__, $prefix);                  } 
	print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Get synthesis of non redundant motifs for $id_bact, ran").": ok\n");     


	if($b_filterout_16SrRNA_res)
	{
	   
	    # --------------------------------------------
	    # filter out motif with 16S RNA sequence
	    my $gb_f = undef;
	    if($id_bact eq $id_bact1){ $gb_f = $gb_f1; }
	    else{                      $gb_f = $gb_f2; }

	    my $dir_of_fasta_for_motifs = $sigffrid_res_dir.'SIGffRid_results/MOTIF_UPSTR_treshmulti'.$MM_order.'_'.$id_bact.'_8/';
	    -e $dir_of_fasta_for_motifs or warn "$prog_tag [WARN] $dir_of_fasta_for_motifs does not exist, mistake in directory name for motifs?, line ".__LINE__."\n";
	    
	    # if no filter is required, file name does not change, otherwise
	    # a new file with name adding _RBSoufiltered before .txt extension 
	    # is created
	    $b_filterout_16SrRNA_res and $grep_sortie = filterout_16SrRNA_res(
		$grep_sortie,
		$gb_f,
		$id_bact,		
		$dir_of_fasta_for_motifs
		);
	    print "$prog_tag $grep_sortie file created\n";
	    print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Filtering of motifs related to 16S sRNA for $id_bact, ran").": ok\n");     
	    # --------------------------------------------	    

	}
	
	# part to sort motifs according to khi2score (usefull??)
	$cmd = "perl ${SIGffRID_scripts_dir}Fsort_by_khi2score.pl -grep_f $grep_sortie";
	$prefix = 'Fsort_by_khi2score_'.short_name($file_with_gene_ids);
	print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
	if($b_sge)     {  push @pids, run_capture::run_capture_sge_prefix($cmd, __LINE__, $prefix);              }
	elsif($b_slurm){  push @pids, run_capture::run_capture_slurm_prefix($cmd, __LINE__, $prefix, $ram_mbgd); }
	else           {  push @pids, run_capture::run_capture_prefix($cmd, __LINE__, $prefix);                  } 
	print(sprintf("%-${display_width_before_column_ok}s", "$prog_tag Motifs sorted by Khi2score for $id_bact, launched").": wait\n");     
    }
    wait_processes_finished(\@pids);
    print "$prog_tag Fsort_by_khi2score part ran.\n$prog_tag All processes finished\n";
    # --------------------------------------------------------

   
    # get end date, something like "Thu Oct 13 04:54:34 1994"
    $end_date = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print "$prog_tag end:$end_date\n";    
}	
	
	
	
