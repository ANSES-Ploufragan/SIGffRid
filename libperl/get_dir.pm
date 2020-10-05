#!/usr/bin/perl
use strict;
use warnings;
$|++;
use Cwd qw(abs_path);

=head1  NAME

get_dir.pm

=head1 DESCRIPTION

Library to get directory (path) name from a file name.

=head1 ARGUMENTS

=over

=item filename <s>        path to deduce directory

=back

example:
get_dir("/home/touzain/dir/file.txt");

=cut

# **********************************************************************

my $prog_tag = '[get_dir.pm]';
my $b_test   = 0; # ok


# ************************************************************************************************
# to extract directory name from file path
# ************************************************************************************************
sub get_dir($)
{
	my($full_path) = @_;
	
	# print "$prog_tag rel path:$full_path becomes ";
	-e $full_path and $full_path = abs_path($full_path);
	$full_path =~ s/^(.+\/).*?$/$1/;
	return $full_path;
}

if($b_test)
{
	my $filename = '/media/sdc_res/raw_sequencage/Mi-seq_Hi-seq/20130712_Miseq_CovetLab_Campilobacter_Integration/assemblies/abyss/CamH0232/CamH0232_S16_L001_R1_001.soft_clean.trimmomatic_trimmed.and_others-unitigs.realign_mapping_on_CamH0232_S16_L001_R1_001.soft_clean.trimmomatic_trimmed.and_others-unitigs_bowtie2-output_f.non-deterministic.bam';
	my $dir      = get_dir($filename);
	print "$prog_tag [Test] dir:$dir extracted from filename:$filename, line ".__LINE__."\n";
}
1;
