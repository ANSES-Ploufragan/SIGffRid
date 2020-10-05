#!/usr/bin/perl
use strict;
use warnings;

package run_capture;
use IPC::System::Simple qw(capture EXIT_ANY);
use Try::Tiny;
use Env qw(ANSES_USER_MAIL SIGFFRID_DIR ANSES_Q_NAME ANSES_SLURM_Q_NAME SLURM_THREADS_A_CORE SLURM_MAX_RAM);

# Copyright or © or Copr. "ANSES Ploufragan, GVB unit  contributor(s) :
# Fabrice Touzain" (2016/01/15)
# 
# fabrice.touzain@anses.fr
# 
# This software is a computer program whose purpose is to [describe
# functionalities and technical features of your software].
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.


=head1  NAME

run_capture.pm

=head1 DESCRIPTION

Library to run command lines from perl using capture tool and sge (multitrheads).

=head1 SUBPROGRAMS

=over

=item run_capture($$;@)

	to run a command line (does not use sge)

=item run_capture_passwd($$;@)

	to run a command line (does not use sge) but
	does not display error to avoid display passwd.
												   
=item run_capture_get_res($$$;@)

	to run a command line and get the result (does
	not use sge)

=item run_capture_sge_prefix($$$;@)

	to run a command line using sge (multithread). Does not wait
	wether job is finished. Third argument allows to give a
	string inserted near the beginning of the script name.

=item run_capture_slurm_prefix($$$;@)

	to run a command line using slurm (multithread). Does not wait
	wether job is finished. Third argument allows to give a
	string inserted near the beginning of the script name.

=item run_capture_sge_prefix_force_exit_status($$$;@)

	to run a command line using sge (multithread). Does not wait
	wether job is finished. Third argument allows to give a
	string inserted near the beginning of the script name.
	Dedicated to command with normal exit status != 0.

=item run_capture_sge($$;@)

	to run a command line using sge (multithread). Does not wait
    wether job is finished.

=item is_run_capture_sge_finished($$)

	to test wether a sge process (using its pid) is finished

=item is_run_capture_slurm_finished($$)

	to test wether a slurm process (using its pid) is finished: MUST BE TESTED

=item is_run_capture_finished($$)

	to test wether a process (using its pid) is finished: MUST BE TESTED
        (job is out of any cluster scheduler)

=item run_capture_get_res_sge_prefix_wait_finished($$$$;@)

	to run command line using sge and get the result,
    wait until job is finished. Third argument allows to give a
    string inserted near the beginning of the script name.
   
=item run_capture_get_res_sge_wait_finished($$$;@)

	to run command line using sge and get the result,
	wait until job is finished.

=item run_capture_sge_prefix_wait_finished($$$;@)

	to run command line using sge (multithread). Wait until job is
	finished. Third argument allows to give a string inserted near
	the beginning of the script name.

=item run_capture_prefix_wait_finished($$$;@)

	to run command line. Wait until job is
	finished. Third argument allows to give a string inserted near
	the beginning of the script name.

=item run_capture_sge_prefix_wait_finished_force_exit_status($$$;@)

	to run command line using sge (multithread). Wait until job is
	finished. Third argument allows to give a string inserted near
	the beginning of the script name.
	Dedicated to command with normal exit status != 0.

=item run_capture_sge_wait_finished($$;@)

	to run command line using sge and wait until job is finished.

=item run_capture_force_ncbi($$;@)

	to run ncbi command line which can fail: try 10 times.

=back

=cut

# **********************************************************************

my $prog_tag                                     = '[run_capture.pm]';
my $SGE_FINISHED_TIME_REQUEST                    = 1; # time between requests (seconds), 5 until 20170511
my $SLURM_FINISHED_TIME_REQUEST                  = 1; # time between requests (seconds), 5 until 20200204
my $FINISHED_TIME_REQUEST                        = 1; # time between requests (seconds), 5 until 20190409
my $b_test_run_capture_sge                       = 0; # ok
my $b_test_run_capture_prefix                    = 1; # 
my $b_test_run_capture_get_res_sge_wait_finished = 0; # ok
# test if sge works for both classic perl call and snakemake call of cluster
my $b_test_run_capture_sge_prefix_wait_finished  = 0; # ok 2018 10 16
my $b_test_run_capture_slurm_prefix_wait_finished= 0; # ok 2020 02 10
my $b_test_run_capture_prefix_wait_finished      = 0; # ok 2019 04 15
# test also good interpretation of cluster command passed via snakemake
my $b_test_run_capture_sge_prefix_wait_finished_force_exit_status = 0; # ok 2018 10 16

# my $ANSES_USER_MAIL                              = $ANSES_USER_MAIL;
if($ANSES_USER_MAIL !~ /^.+\@.+$/)
{
	die "$prog_tag [Error] Invalid ANSES_USER_MAIL environment variable to mail cluster info, please check, line ".__LINE__."\n";
}

# directory for cluster errors and outputs
my $op_dir   = $SIGFFRID_DIR.'/cluster_errors/';
my $err_dir  = $op_dir;
-e $op_dir or die "$prog_tag [Error] $op_dir output directory for SGE does not exist, please check, line ".__LINE__."\n";

my $STD_SCRIPT_HEAD_SGE = "
#!/bin/bash

# Les commentraires qui commencent par '#\$' sont
# interpretes par SGE comme des options en ligne

# Shell a utiliser pour l'execution du job
#\$ -S /bin/bash

# Export de toutes les variables d'environnement
#\$ -V

# Utilisateur a avertir, fabrice.touzain\@anses.fr
#\$ -M $ANSES_USER_MAIL

# Avertir au debut (b)egin, a la fin (e)nd, a echec (a)bort et
# a la suspension (s)uspend d'un job
#\$ -m ea

# Sortie standard
# Vous pouvez utiliser '-j y' pour ajouter stderr avec stdout
#\$ -o $op_dir

# Sortie d'erreur (ne pas utiliser cette option avec '-j y')
#\$ -e $err_dir

# Lance la commande depuis le repertoire ou est lance le script
#\$ -cwd
# Exemple pour l'utilisation de rmes
    ";

# slurm options
my $SLURM_ntasks_per_node = $SLURM_THREADS_A_CORE - 1;
# my $SLURM_MAX_RAM         = 28672; # ENV VAR
my $SLURM_MAX_TIME_DAYS   = 10;
my $SLURM_MAX_TIME_HOURS  = 00;
# ---------------------------------------------------------------------------------

=head2  NAME

untaint_cmd($)

=head3 DESCRIPTION

Check that cmd string has only safe characters not to trigger dangerous commands.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)

=back

=head3 RETURN

=over

=item safe string to use as a command.

=back

=cut

sub untaint_cmd($)
{
	my($string) = @_;
	
	if ($string =~ /^([\-\@\w\d_\[\]'\/\s\*\.\:\%\>\^]+)$/) {
	    $string = $1;                     # $string now untainted.
	}
	else {
		my $forbidden = $string;
		$forbidden =~ s/[\-\@\w\d_\[\]'\/\s\*\.\:\%\>\^]+//g;
	    die "$prog_tag [untaint_cmd][Error] Bad data ($forbidden) in '$string', line ".__LINE__."\n";        # Log this somewhere.
	}
	return $string;
}
# ---------------------------------------------------------------------------------

=head2  NAME

run_capture($$;@)

=head3 DESCRIPTION

to run linux command through capture perl lib when we dont want get output
(does not use sge).

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=cut

sub run_capture($$;@)
{
	my($cmd,$line,@error_allowed) = @_;

	# $cmd = untaint_cmd($cmd);
	
	my $msg = '';
	my( $success, $full_buf, $stdout_buf, $stderr_buf );
	
	# run here
	try
	{
		if(@error_allowed)
		{
		    if($error_allowed[0] eq 'EXIT_ANY')
			{
			    $msg = capture(EXIT_ANY,$cmd);
			}
			else
			{
			    $msg = capture([@error_allowed],$cmd);
			}
		}
		else
		{
			$msg = capture($cmd);
		}
	}
	catch
	{
		# die "$prog_tag [Error] cmd:$cmd failed, msg:$msg, line $line, error:".join('',@$stderr_buf)."\n";
		die "$prog_tag [Error] cmd:$cmd failed, msg:$msg, line $line, error:$@\n";
	};
	print "$msg\n";
}
# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_passwd($$;@)

=head3 DESCRIPTION

to run linux command through capture perl lib when we dont want get output
(does not use sge) and do not want error message with string written by user.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture_passwd($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=cut

sub run_capture_passwd($$;@)
{
	my($cmd,$line,@error_allowed) = @_;

	# $cmd = untaint_cmd($cmd);
	
	my $msg = '';
	# run here
	try
	{
		if(@error_allowed)
		{
		        if($error_allowed[0] eq 'EXIT_ANY')
			  {
			    $msg = capture(EXIT_ANY,$cmd);
			  }
			else
			  {
			    $msg = capture([@error_allowed],$cmd);
			  }
		}
		else
		{
			$msg = capture($cmd);
		}
	}
	catch
	{
		die "$prog_tag [Error] cmd failed, probably bad passwd\n";
	};
	print "$msg\n";
}
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_get_res($$$;@)

=head3 DESCRIPTION

to run linux command through capture perl lib when we want get output
(does not use sge).

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item r_msg: ref to an array where output will be stored

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=cut

sub run_capture_get_res($$$;@)
{
	my($cmd,$line,$r_msg,@error_allowed) = @_;

	# $cmd = untaint_cmd($cmd);
	
	# run here
	try
	{
		if(@error_allowed)
		{
		  		
		    if($error_allowed[0] eq 'EXIT_ANY')
		    {
			# print "$prog_tag [run_capture_get_res] try capture with EXIT_ANY, line ".__LINE__."\n";
			@$r_msg = capture(EXIT_ANY,$cmd);
		    }
		    else
		    {
			# print "$prog_tag [run_capture_get_res] try capture with error_allowed array, line ".__LINE__."\n";
			@$r_msg = capture([@error_allowed],$cmd);
		    }
		}
		else
		{
		    # print "$prog_tag [run_capture_get_res] try single capture, line ".__LINE__."\n";
		    @$r_msg = capture($cmd);
		}
	}
	catch
	{
	    die "$prog_tag [run_capture_get_res][Error] cmd:$cmd failed, line $line, curr_line ".__LINE__."\n".join("\n",'msg:',@$r_msg,",error:$@",'');
	};
}
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_sge_prefix($$$;@)

=head3 DESCRIPTION

FUNCTION: to run linux command through capture perl lib using SGE (cluster).

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script file name

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=head3 RETURN

the pid of the cluster process

=head3 DEPENDS ON

run_capture_get_res($$$;@)

=cut

sub run_capture_sge_prefix($$$;@)
{
	# use Cwd;
	use File::Temp qw/ tempfile /;
	
	my($cmd,$line,$prefix,@error_allowed) = @_;

	my $pid = undef;
	my @res = ();

	## sge options are stored in this file
	#my $sge_option_file = "run_capture_sge_option_file.txt";
	#-e $sge_option_file or die "$prog_tag [Error] [run_capture_sge] $sge_option_file file does not exist, line ".__LINE__.", run_capture_sge called from line $line\n";

	# sh script file submitted to sge, create as a tmp file
	my ($fh, $sge_sh_scripts_f) = tempfile('run_capture_sge_'.$prefix.'_XXXXX',
											DIR    => $op_dir,
											UNLINK => 0,
											SUFFIX => '.sh');
	binmode( $fh, ":utf8" );

	my $sge_sh_scripts = $STD_SCRIPT_HEAD_SGE."$cmd\n";
	
	# print sh scripts to run in a file
	print $fh $sge_sh_scripts;
	close $fh;
	
	# puts appropriate executive rights to files (run...format)
	my $rights = 0755;
	chmod($rights,$sge_sh_scripts_f) or die "$prog_tag [Error] [run_capture_sge_prefix] Cannot chmod $rights for $sge_sh_scripts_f file:$!, line ".__LINE__."\n";
	
	# add sge submission command with option parameter file and the command asked by user
	my $sge_cmd = 'qsub -q '.$ANSES_Q_NAME." $sge_sh_scripts_f";

	# run sge command in shell
	print "$prog_tag [run_capture_sge_prefix] sge_cmd:$sge_cmd, line ".__LINE__."\n";
	run_capture_get_res($sge_cmd,$line,\@res, @error_allowed);
	print "$prog_tag [run_capture_sge_prefix] ran\n";

	# print "res:".join(',',@res)."\n";

	foreach(@res)
	{
		/^Your job (\d+)/ and do
		{	
			$pid = $1;
			unlink $sge_sh_scripts_f;
			last;
		};
	}
	if(not defined $pid)
	{
		die join('', "$prog_tag [Error] [run_capture_sge_prefix] We got something else than pid:\n",
			@res, ", line ".__LINE__.", run_capture_sge called from line $line\n");
	}
	return($pid);
}
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_slurm_prefix($$$;@)

=head3 DESCRIPTION

FUNCTION: to run linux command through capture perl lib using SLURM (cluster).

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script file name

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=head3 RETURN

the pid of the cluster process

=head3 DEPENDS ON

run_capture_get_res($$$;@)

=cut

sub run_capture_slurm_prefix($$$$;@)
{
	# use Cwd;
	use File::Temp qw/ tempfile /;
	
	my($cmd,$line,$prefix,$job_mem, @error_allowed) = @_;

	my $pid = undef;
	my @res = ();

	## slurm options are stored in this file
	#my $slurm_option_file = "run_capture_slurm_option_file.txt";
	#-e $slurm_option_file or die "$prog_tag [Error] [run_capture_slurm] $slurm_option_file file does not exist, line ".__LINE__.", run_capture_slurm called from line $line\n";


	my $STD_SCRIPT_HEAD_SLURM = "#!/bin/bash
#
#SBATCH --mail-type=END
#SBATCH --mail-user=$ANSES_USER_MAIL
#SBATCH -p $ANSES_SLURM_Q_NAME # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --tasks-per-node 1 # $SLURM_ntasks_per_node
#SBATCH --mem $job_mem # $SLURM_MAX_RAM # memory pool for all cores
#SBATCH --mpi=pmix_v2
#SBATCH -t $SLURM_MAX_TIME_DAYS-$SLURM_MAX_TIME_HOURS:00 # time (D-HH:MM)
#SBATCH -o ${op_dir}${prefix}_slurm.%N.%j.out # STDOUT
#SBATCH -e ${op_dir}${prefix}_slurm.%N.%j.err # STDERR
    ";

	# sh script file submitted to slurm, create as a tmp file
	my ($fh, $slurm_sh_scripts_f) = tempfile('run_capture_slurm_'.$prefix.'_XXXXX',
											DIR    => $op_dir,
											UNLINK => 0,
											SUFFIX => '.sh');
	binmode( $fh, ":utf8" );

	my $slurm_sh_scripts = $STD_SCRIPT_HEAD_SLURM."$cmd\n";
	
	# print sh scripts to run in a file
	print $fh $slurm_sh_scripts;
	close $fh;
	
	# puts appropriate executive rights to files (run...format)
	my $rights = 0755;
	chmod($rights,$slurm_sh_scripts_f) or die "$prog_tag [Error] [run_capture_slurm_prefix] Cannot chmod $rights for $slurm_sh_scripts_f file:$!, line ".__LINE__."\n";
	
	# add slurm submission command with option parameter file and the command asked by user
	my $slurm_cmd = "srun --mpi=pmi2 $slurm_sh_scripts_f";

	# run slurm command in shell
	print "$prog_tag [run_capture_slurm_prefix] slurm_cmd:$slurm_cmd, line ".__LINE__."\n";
	run_capture_get_res($slurm_cmd,$line,\@res, @error_allowed);
	print "$prog_tag [run_capture_slurm_prefix] ran\n";

	# print "res:".join(',',@res)."\n";

	foreach(@res)
	{
		/^Submitted batch job (\d+)/ and do
		{	
			$pid = $1;
			unlink $slurm_sh_scripts_f;
			last;
		};
	}
	if(not defined $pid)
	{
		die join('', "$prog_tag [Error] [run_capture_slurm_prefix] We got something else than pid:\n",
			@res, ", line ".__LINE__.", run_capture_slurm called from line $line\n");
	}
	return($pid);
}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_prefix($$$;@)

=head3 DESCRIPTION

FUNCTION: to run linux command through capture perl lib using SGE (cluster).

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script file name

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=head3 RETURN

the pid of the cluster process

=head3 DEPENDS ON

run_capture_get_res($$$;@)

=cut

sub run_capture_prefix($$$;@)
{
    use Proc::Daemon;
    # use Cwd 'abs_path';
    # use IPC::Run3 qw(run3);
    use File::Temp qw/ tempfile /;
    
    my($cmd,$line,$prefix,@error_allowed) = @_;
    
    my $pid = undef;
    my @res = ();
    
    # sh script file submitted to bash, create as a tmp file
    my ($fh, $sh_scripts_f) = tempfile('run_capture_'.$prefix.'_XXXXX',
				       DIR    => $op_dir,
				       UNLINK => 0,
				       SUFFIX => '.sh');
    binmode( $fh, ":utf8" );
    
    my $error_f = $sh_scripts_f;
    my $out_f = $sh_scripts_f;	
    $error_f =~ s/\.sh$/.e/;
    $out_f =~ s/\.sh$/.o/;	
    
    my $BASH_SHEBANG = "#!/bin/bash";
    my $sh_scripts = "$BASH_SHEBANG\n$cmd\n";
    
    # print sh scripts to run in a file
    print $fh $sh_scripts;
    close $fh;
    
    # puts appropriate executive rights to files (run...format)
    my $rights = 0755;
    chmod($rights,$sh_scripts_f) or die "$prog_tag [Error] [run_capture_prefix] Cannot chmod $rights for $sh_scripts_f file:$!, line ".__LINE__."\n";
    
    # add sge submission command with option parameter file and the command asked by user
    my $sh_cmd = "$sh_scripts_f 1> $out_f 2> $error_f";
    
    print "$prog_tag [run_capture_prefix] sh_cmd:$sh_cmd, line ".__LINE__."\n";
    
    # open my $fho, ">", $out_f;
    # open my $fhe, ">", $error_f;	
    # run3([$sh_scripts_f], $error_f, $out_f);
    # close $fho;
    # close $fhe;

    my $daemon = Proc::Daemon->new();
    $pid = $daemon->Init( {
    	work_dir     => $op_dir,
    	exec_command => [ $sh_scripts_f ],
    	child_STDOUT => $out_f,
    	child_STDERR => $error_f
    			  } );
    
    print "$prog_tag [run_capture_prefix] launched with pid:$pid, $out_f and $error_f files created\n";
    
    return($pid);
}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_sge_prefix_force_exit_status($$$;@)

=head3 DESCRIPTION

FUNCTION: to run linux command through capture perl lib using SGE (cluster).

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script file name

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=head3 RETURN

the pid of the cluster process

=head3 DEPENDS ON

run_capture_get_res($$$;@)

=cut

sub run_capture_sge_prefix_force_exit_status($$$;@)
{
	# use Cwd;
	use File::Temp qw/ tempfile /;
	
	my($cmd,$line,$prefix,@error_allowed) = @_;

	my $pid = undef;
	my @res = ();

	## sge options are stored in this file
	#my $sge_option_file = "run_capture_sge_option_file.txt";
	#-e $sge_option_file or die "$prog_tag [Error] [run_capture_sge] $sge_option_file file does not exist, line ".__LINE__.", run_capture_sge called from line $line\n";


	# sh script file submitted to sge, create as a tmp file
	my ($fh, $sge_sh_scripts_f) = tempfile('run_capture_sge_'.$prefix.'_XXXXX',
											DIR    => $op_dir,
											UNLINK => 0,
											SUFFIX => '.sh');
	binmode( $fh, ":utf8" );

	my $sge_sh_scripts = $STD_SCRIPT_HEAD_SGE."$cmd\n";
	
	# print sh scripts to run in a file
	print $fh $sge_sh_scripts;
	close $fh;
	
	# puts appropriate executive rights to files (run...format)
	my $rights = 0755;
	chmod($rights,$sge_sh_scripts_f) or die "$prog_tag [Error] [run_capture_sge_prefix] Cannot chmod $rights for $sge_sh_scripts_f file:$!, line ".__LINE__."\n";
	
	# add sge submission command with option parameter file and the command asked by user
	my $sge_cmd = 'qsub -sync y -q '.$ANSES_Q_NAME." $sge_sh_scripts_f";

	# run sge command in shell
	print "$prog_tag [run_capture_sge_prefix] sge_cmd:$sge_cmd, line ".__LINE__."\n";
	run_capture_get_res($sge_cmd,$line,\@res, ('EXIT_ANY'));
	print "$prog_tag [run_capture_sge_prefix] ran\n";

	# print "res:".join(',',@res)."\n";

	foreach(@res)
	{
		/^Your job (\d+)/ and do
		{	
			$pid = $1;
			unlink $sge_sh_scripts_f;
			last;
		};
	}
	if(not defined $pid)
	{
		die join('', "$prog_tag [Error] [run_capture_sge_prefix] We got something else than pid:\n",
			@res, ", line ".__LINE__.", run_capture_sge called from line $line\n");
	}
	#if($res[0] =~ /^Your job (\d+)/)
	#{
		#$pid = $1;
		#unlink $sge_sh_scripts_f;
	#}
	#else
	#{
		#die "$prog_tag [Error] [run_capture_sge_prefix] We got something else than pid ($res[0]), line ".__LINE__.", run_capture_sge called from line $line\n";
	#}
	return($pid);
}
# ---------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_sge($$;@)

=head3 DESCRIPTION

Run capture_sge_prefix with an empty prefix.
Made to ensure compatibility.

=cut

sub run_capture_sge($$;@)
{
	my @options = @_;
	if(scalar(@options)>2)
	{
		run_capture_sge_prefix($options[0], $options[1],'unkowntask_',@options[2..$#options]);
	}
	else
	{
		run_capture_sge_prefix($options[0], $options[1],'unkowntask_');
	}
}
# ---------------------------------------------------------------------------------

=head2  NAME

is_run_capture_sge_finished($$)

=head3 DESCRIPTION

FUNCTION: test wether the process with a given pid is finished on the cluster

=head3 ARGUMENTS

=over

=item pid:  number of the process to check (integer)

=item line: the line (location) where this sub is called (__LINE__)

=back

=head3 GLOBAL VAR

SGE_FINISHED_TIME_REQUEST:5 seconds (times between two requests to qstat)

=head3 DEPENDS ON

run_capture_get_res($$$;@)

=head3 RETURN

A boolean: 1 or 0

=cut

sub is_run_capture_sge_finished($$)
{
	my($pid,$line) = @_;

	defined $pid or return 1;
	
	my $cmd        = "qstat -j $pid";
	my @res        = ();
	my $b_finished = 0;

	# run qstat to know the state of the pid
	run_capture_get_res($cmd,$line,\@res,'EXIT_ANY');

	# means the job is already finished
	# despite what is written on screen 'Following jobs do not exist' than pid on following line, nothing is returned
	if($#res == -1)
	{
		$b_finished = 1;
	}
	else
	{
		foreach(@res)
		{
			/Following jobs do not exist/ and do
			{
				$b_finished = 1;
				last;
			};
		}
		## means something else
		#foreach($res[0])
		#{
			## cases to add if you want a special messages for some errors
			#print "$prog_tag [Warn][run_capture_sge] Cluster message unknown:\n".join('',@res).", line ".__LINE__.", run_capture_sge called line $line\n";
			#last;
		#}

		# job is not finished but in case of cluster error, we can personalize msg in commented loop above 
		# $b_finished = 0;
	}
	
	return($b_finished);
}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

is_run_capture_slurm_finished($$)

=head3 DESCRIPTION

FUNCTION: test wether the process with a given pid is finished on the cluster

=head3 ARGUMENTS

=over

=item pid:  number of the process to check (integer)

=item line: the line (location) where this sub is called (__LINE__)

=back

=head3 GLOBAL VAR

SLURM_FINISHED_TIME_REQUEST:5 seconds (times between two requests to qstat)

=head3 DEPENDS ON

run_capture_get_res($$$;@)

=head3 RETURN

A boolean: 1 or 0

=cut

sub is_run_capture_slurm_finished($$)
{
	my($pid,$line) = @_;

	defined $pid or return 1;
	
	my $cmd        = "squeue -j$pid";
	my @res        = ();
	my $b_finished = 0;

	# run qstat to know the state of the pid
	run_capture_get_res($cmd,$line,\@res,'EXIT_ANY');

	# means the job is already finished
	# always 1 line
	# JOBID PARTITION ...
	if($#res == 0)
	{
		$b_finished = 1;
	}
	else
	{
		foreach(@res)
		{
			/^\s*$pid[\s\t]+/ and do
			{
				$b_finished = 0;
				last;
			};
		}
		## means something else
		#foreach($res[0])
		#{
			## cases to add if you want a special messages for some errors
			#print "$prog_tag [Warn][run_capture_slurm] Cluster message unknown:\n".join('',@res).", line ".__LINE__.", run_capture_slurm called line $line\n";
			#last;
		#}

		# job is not finished but in case of cluster error, we can personalize msg in commented loop above 
		# $b_finished = 0;
	}
	
	return($b_finished);
}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

is_run_capture_finished($$)

=head3 DESCRIPTION

FUNCTION: test wether the process with a given pid is finished on the cluster

=head3 ARGUMENTS

=over

=item pid:  number of the process to check (integer)

=item line: the line (location) where this sub is called (__LINE__)

=back

=head3 GLOBAL VAR

FINISHED_TIME_REQUEST:5 seconds (times between two requests to qstat)

=head3 DEPENDS ON

run_capture_get_res($$$;@)

=head3 RETURN

A boolean: 1 or 0

=cut

sub is_run_capture_finished($$)
{
	my($pid,$line) = @_;

	defined $pid or return 1;
	
	my $cmd        = "ps -p $pid";
	my @res        = ();
	my $b_finished = 1;

	# run qstat to know the state of the pid
	run_capture_get_res($cmd,$line,\@res,'EXIT_ANY');

	# means the job is already finished
	# you always have at less 1 line:
	#   PID TTY          TIME CMD
	if($#res == 0)
	{
		$b_finished = 1;
	}
	else
	{
	    # typical output:
	    #   PID TTY          TIME CMD
	    #24871 pts/21   00:00:00 telnet 
	    foreach(@res)
	    {
		/^\s*$pid / and do
		{
		    $b_finished = 0;
		    last;
		};
	    }
	    ## means something else
	    #foreach($res[0])
	    #{
	    ## cases to add if you want a special messages for some errors
	    #print "$prog_tag [Warn][run_capture_slurm] Cluster message unknown:\n".join('',@res).", line ".__LINE__.", run_capture_slurm called line $line\n";
	    #last;
	    #}
	    
	    # job is not finished but in case of cluster error, we can personalize msg in commented loop above 
	    # $b_finished = 0;
	}
	
	return($b_finished);
}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

are_half_of_these_pid_finished($$$)

=head3 DESCRIPTION

FUNCTION: test wether half of the processes with given pids are finished on the cluster.

=head3 ARGUMENTS

=over

=item r_pid: ref to an array of pids of the process to check (integer)

=item line: the line (location) where this sub is called (__LINE__)

=back

=head3 GLOBAL VAR

SGE_FINISHED_TIME_REQUEST:5 seconds (times between two requests to qstat)

=head3 DEPENDS ON

=over

=item run_capture_get_res($$$;@)

=item is_run_capture_sge_finished($$)

=back

=head3 RETURN

A boolean: 1 or 0

=cut

sub are_unfinished_pid_less_than_half_avail_threads($$$)
{
	my($r_pid,$nb_threads,$line) = @_;
	my $nb_finished          = 0;

	# this variable tells than max 16/max_ratio_of_threads assemblies
	# main processes can be running simultaneously (4)
	my $max_ratio_of_threads = 4;

	my $nb_pids          = scalar(@$r_pid);
	$nb_pids or return 1;

	use POSIX qw(ceil);
	my $nb_threads_overn = ceil($nb_threads /$max_ratio_of_threads);

	my @pids_to_keep = ();
	foreach my $pid(@$r_pid)
	{
		my $cmd = "qstat -j $pid";
		my @res = ();
		
	
		# run qstat to know the state of the pid
		run_capture_get_res($cmd,$line,\@res,'EXIT_ANY');
	
		# means the job is already finished
		# despite what is written on screen 'Following jobs do not exist' than pid on following line, nothing is returned
		if($#res == -1)
		{
			$nb_finished++;
			# next;
		}
		else
		{
			push @pids_to_keep, $pid;
			
			## means something else
			#foreach($res[0])
			#{
				## cases to add if you want a special messages for some errors
				#print "$prog_tag [Warn][run_capture_sge] Cluster message unknown:\n".join('',@res).", line ".__LINE__.", run_capture_sge called line $line\n";
				#last;
			#}
	
			# job is not finished but in case of cluster error, we can personalize msg in commented loop above 
			# $b_finished = 0;
		}
	}

	# we keep only not finished jobs to avoid to list all finished processes all the time
	@$r_pid = @pids_to_keep;
		
	return(($nb_pids - $nb_finished) < $nb_threads_overn);
}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_sge_prefix_wait_halfunfinishedpids($$$$$;@)

=head3 DESCRIPTION

FUNCTION: to run linux command through capture perl lib using SGE (cluster).
Before launching the command, wait that number of pids still running in the
argument list is lower than half of the available theads.


=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script file name

=item reference to an array of ran pids

=item number of threads

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=head3 RETURN

the pid of the cluster process

=head3 DEPENDS ON

(run_capture_get_res($$$;@))
run_capture_sge_prefix($$$;@)
are_unfinished_pid_less_than_half_avail_threads($$$)

=cut

sub run_capture_sge_prefix_wait_halfunfinishedpids($$$$$;@)
{
	my($cmd,$line,$prefix,$r_assemby_pids_launched,$nb_threads,@error_allowed) = @_;
	
	# ------------------------------------------------------------------
	# wait that number of pids still running in the
	# argument list is lower than half of the available theads
	# ------------------------------------------------------------------
	while(not are_unfinished_pid_less_than_half_avail_threads($r_assemby_pids_launched,$nb_threads,$line))
	{
		sleep $SGE_FINISHED_TIME_REQUEST;
	}
	# ------------------------------------------------------------------

	# run process on the cluster
	run_capture_sge_prefix($cmd,$line,$prefix,@error_allowed);
}	
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

are_unfinished_pid_less_than_nb_max_pids($$$)

=head3 DESCRIPTION

FUNCTION: test weather number of the processes with given pids is less than
the number of allowed threads for the global work launch requested on the
cluster.

=head3 ARGUMENTS

=over

=item r_pid: ref to an array of pids of the process to check (integer)

=item line: the line (location) where this sub is called (__LINE__)

=back

=head3 GLOBAL VAR

SGE_FINISHED_TIME_REQUEST:5 seconds (times between two requests to qstat)

=head3 DEPENDS ON

=over

=item run_capture_get_res($$$;@)

=item is_run_capture_sge_finished($$)

=back

=head3 RETURN

A boolean: 1 or 0

=cut

my $b_test_are_unfinished_pid_less_than_nb_max_pids = 0; # OK 
sub are_unfinished_pid_less_than_nb_max_pids($$$)
{
	my($r_pid,$nb_max_pids,$line) = @_;
	my $nb_finished          = 0;

	my $nb_pids          = scalar(@$r_pid);
	
	if($nb_pids == 0)
	{
		$b_test_are_unfinished_pid_less_than_nb_max_pids and
			print "$prog_tag [are_unfinished_pid_less_than_nb_max_pids] nb pids remaining: $nb_pids (".join(':',@$r_pid)."), return 1, line ".__LINE__."\n";
		return 1;
	}
	
	my @pids_to_keep = ();
	#foreach my $pid(@$r_pid)
	#{
		#my $cmd = "qstat -j $pid";
		#my @res = ();
		
	
		## run qstat to know the state of the pid
		#run_capture_get_res($cmd,$line,\@res,'EXIT_ANY');
	
		## means the job is not finished
		## despite what is written on screen 'Following jobs do not exist'
		## than pid on following line, nothing is returned when a job is
		## finished
		#if($#res != -1)
		#{
			#push @pids_to_keep, $pid;
		#}
	#}
	
	# keeps only IDs of current running jobs, no need loop
	my $cmd = "qstat -j ".join(',', @$r_pid)." | grep job_number | cut -d ' ' -f 2- | tr -d ' ' ";
	my @res = ();
	
	# run qstat to know the state of the pid
	run_capture_get_res($cmd,$line,\@res,'EXIT_ANY');

	# means the job is not finished
	# despite what is written on screen 'Following jobs do not exist'
	# than pid on following line, nothing is returned when a job is
	# finished
	if($#res != -1)
	{
		chomp(@res);
		push @pids_to_keep, @res;
	}

	# we keep only not finished jobs to avoid to list all finished processes all the time
	@$r_pid = @pids_to_keep;
	$nb_pids = scalar(@pids_to_keep);
	$b_test_are_unfinished_pid_less_than_nb_max_pids and
		print "$prog_tag [are_unfinished_pid_less_than_nb_max_pids] nb pids remaining: $nb_pids (".join(':',@$r_pid)."), return ".($nb_pids < $nb_max_pids).",line ".__LINE__."\n";

	if($nb_pids < $nb_max_pids)
	{
		$b_test_are_unfinished_pid_less_than_nb_max_pids and
			print "$prog_tag [are_unfinished_pid_less_than_nb_max_pids] RETURN 1, line ".__LINE__."\n";
		return 1;
	}
	else
	{
		$b_test_are_unfinished_pid_less_than_nb_max_pids and
			print "$prog_tag [are_unfinished_pid_less_than_nb_max_pids] RETURN 0, line ".__LINE__."\n";
		return 0;
	}
	# return($nb_pids < $nb_max_pids);
}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_sge_prefix_wait_maxnbpids($$$$$;@)

=head3 DESCRIPTION

FUNCTION: to run linux command through capture perl lib using SGE (cluster).
Before launching the command, wait that number of pids still running in the
argument list is lower than max number of theads allowed by user.


=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string) run_capture($$;@)                            to run a command line
												   (does not use sge)

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script file name

=item reference to an array of ran pids

=item max number of allowed threads simultaneously

=item errors_allowed: list of errors allowed (and caught), EXIT_ANY for all.

=back

=head3 RETURN

the pid of the cluster process

=head3 DEPENDS ON

(run_capture_get_res($$$;@))
run_capture_sge_prefix($$$;@)
are_unfinished_pid_less_than_nb_max_pids($$$)

=cut

sub run_capture_sge_prefix_wait_maxnbpids($$$$$;@)
{
	my($cmd,$line,$prefix,$r_assemby_pids_launched,$nb_max_pids,@error_allowed) = @_;
	
	# ------------------------------------------------------------------
	# wait that number of pids still running in the
	# argument list is lower than half of the available theads
	# ------------------------------------------------------------------
	while(not are_unfinished_pid_less_than_nb_max_pids($r_assemby_pids_launched,$nb_max_pids,$line))
	{
		sleep $SGE_FINISHED_TIME_REQUEST;
	}
	# ------------------------------------------------------------------

	# run process on the cluster
	push @$r_assemby_pids_launched, run_capture_sge_prefix($cmd,$line,$prefix,@error_allowed);
}	
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_get_res_sge_prefix_wait_finished($$$$;@)

=head3 DESCRIPTION

To run linux command through capture perl lib using SGE (cluster)
when we want to get results.
Do not give back the hand until process is finished.
Third argument gives a string inserted near the beginning of the process (script) name.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use sge

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: the string inserted near the beginning of the script name ran.

=item r_msg: ref to an array where we store output (got from cluster out
      in fact)
      
=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 DEPENDS ON

=over

=item is_run_capture_sge_finished($$)

=item run_capture_get_res($$$;@)

=back

=cut

sub run_capture_get_res_sge_prefix_wait_finished($$$$;@)
{
	# use Cwd;
	use File::Temp qw/ tempfile /;
	
	my($cmd,$line,$prefix,$r_msg,@error_allowed) = @_;
	-e $SIGFFRID_DIR or die "$prog_tag [run_capture_get_res_sge_wait_finished][Error] SIGFFRID_DIR environment variable not available or wrong ($SIGFFRID_DIR), line ".__LINE__."\n";

	my $pid = undef;
	my @res = ();

	## sge options are stored in this file
	#my $sge_option_file = "run_capture_sge_option_file.txt";
	#-e $sge_option_file or die "$prog_tag [Error] [run_capture_sge] $sge_option_file file does not exist, line ".__LINE__.", run_capture_sge called from line $line\n";


	# sh script file submitted to sge, create as a tmp file
	my ($fh, $sge_sh_scripts_f) = tempfile('run_capture_sge_'.$prefix.'_XXXXX',
											DIR    => $op_dir,
											UNLINK => 0,
											SUFFIX => '.sh');
	binmode( $fh, ":utf8" );

	my $sge_sh_scripts   = $STD_SCRIPT_HEAD_SGE."$cmd\n";
	
	# print sh scripts to run in a file
	print $fh $sge_sh_scripts;
	close $fh;
	
	# puts appropriate executive rights to files (run...format)
	my $rights = 0755;
	chmod($rights,$sge_sh_scripts_f) or die "$prog_tag [run_capture_get_res_sge_prefix_wait_finished][Error] Cannot chmod $rights for $sge_sh_scripts_f file:$!\n";
	
	# add sge submission command with option parameter file and the command asked by user
	my $sge_cmd = "qsub -q $ANSES_Q_NAME -p -1 $sge_sh_scripts_f";

	# run sge command in shell
	print "$prog_tag [run_capture_get_res_sge_wait_finished] sge_cmd:$sge_cmd, line ".__LINE__."\n";
	run_capture_get_res($sge_cmd,$line,\@res, @error_allowed);
	print "$prog_tag [run_capture_get_res_sge_wait_finished] ran\n";

	print "res:".join(',',@res)."\n";

	if($res[$#res] =~ /^Your job (\d+)/)
	{
		$pid = $1;
		unlink $sge_sh_scripts_f;
	}
	else
	{
		die "$prog_tag [run_capture_get_res_sge_wait_finished][Error] We got something else than pid ($res[$#res]), line ".__LINE__.", run_capture_sge called from line $line\n";
	}

	# we wait that cluster job is finnished
	while(not is_run_capture_sge_finished($pid,$line))
	{
		sleep $SGE_FINISHED_TIME_REQUEST;
	}
	
	# now we can get output in the dedicated cluster file

	# get the name of the process output file
	my @cluster_out_f = glob($SIGFFRID_DIR."cluster_errors/*.o$pid");

	my $nb_f = scalar(@cluster_out_f);
	if($nb_f == 1)
	{
		# we get output and store it in the array given as argument
		open(SCOF,'<',$cluster_out_f[0])or die "$prog_tag [Error] Cannot open $cluster_out_f[0] file:$!, line ".__LINE__."\n";
		while(<SCOF>)
		{
			push @$r_msg, $_;
		}
		close SCOF;
	}
	else
	{
		die "$prog_tag [run_capture_get_res_sge_wait_finished][Error] We find $nb_f cluster output for pid:$pid, we might find only one, line ".__LINE__."\n";
	}
}
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_get_res_slurm_prefix_wait_finished($$$$;@)

=head3 DESCRIPTION

To run linux command through capture perl lib using SLURM (cluster)
when we want to get results.
Do not give back the hand until process is finished.
Third argument gives a string inserted near the beginning of the process (script) name.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use slurm

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: the string inserted near the beginning of the script name ran.

=item r_msg: ref to an array where we store output (got from cluster out
      in fact)
      
=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 DEPENDS ON

=over

=item is_run_capture_slurm_finished($$)

=item run_capture_get_res($$$;@)

=back

=cut

sub run_capture_get_res_slurm_prefix_wait_finished($$$$$;@)
{
	# use Cwd;
	use File::Temp qw/ tempfile /;
	
	my($cmd,$line,$prefix,$r_msg,$job_mem,@error_allowed) = @_;
	-e $SIGFFRID_DIR or die "$prog_tag [run_capture_get_res_slurm_wait_finished][Error] SIGFFRID_DIR environment variable not available or wrong ($SIGFFRID_DIR), line ".__LINE__."\n";

	my $pid = undef;
	my @res = ();

	my $STD_SCRIPT_HEAD_SLURM = "#!/bin/bash
#
#SBATCH --mail-type=END
#SBATCH --mail-user=$ANSES_USER_MAIL
#SBATCH -p $ANSES_SLURM_Q_NAME # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mpi=pmix_v2
#SBATCH --tasks-per-node 1 # $SLURM_ntasks_per_node
#SBATCH --mem $job_mem # $SLURM_MAX_RAM # memory pool for all cores
#SBATCH -t $SLURM_MAX_TIME_DAYS-$SLURM_MAX_TIME_HOURS:00 # time (D-HH:MM)
#SBATCH -o ${op_dir}${prefix}_slurm.%N.%j.out # STDOUT
#SBATCH -e ${op_dir}${prefix}_slurm.%N.%j.err # STDERR
    ";


	## slurm options are stored in this file
	#my $slurm_option_file = "run_capture_slurm_option_file.txt";
	#-e $slurm_option_file or die "$prog_tag [Error] [run_capture_slurm] $slurm_option_file file does not exist, line ".__LINE__.", run_capture_slurm called from line $line\n";


	# sh script file submitted to slurm, create as a tmp file
	my ($fh, $slurm_sh_scripts_f) = tempfile('run_capture_slurm_'.$prefix.'_XXXXX',
						 DIR    => $op_dir,
						 UNLINK => 0,
						 SUFFIX => '.sh');
	binmode( $fh, ":utf8" );

	my $slurm_sh_scripts   = $STD_SCRIPT_HEAD_SLURM."$cmd\n";
	
	# print sh scripts to run in a file
	print $fh $slurm_sh_scripts;
	close $fh;
	
	# puts appropriate executive rights to files (run...format)
	my $rights = 0755;
	chmod($rights,$slurm_sh_scripts_f) or die "$prog_tag [run_capture_get_res_slurm_prefix_wait_finished][Error] Cannot chmod $rights for $slurm_sh_scripts_f file:$!\n";
	
	# add slurm submission command with option parameter file and the command asked by user
	my $slurm_cmd = "srun -mpi=pmi2 $slurm_sh_scripts_f";

	# run slurm command in shell
	print "$prog_tag [run_capture_get_res_slurm_wait_finished] slurm_cmd:$slurm_cmd, line ".__LINE__."\n";
	run_capture_get_res($slurm_cmd,$line,\@res, @error_allowed);
	print "$prog_tag [run_capture_get_res_slurm_wait_finished] ran\n";

	print "res:".join(',',@res)."\n";

	if($res[$#res] =~ /^Submitted batch job (\d+)/)
	{
		$pid = $1;
		unlink $slurm_sh_scripts_f;
	}
	else
	{
		die "$prog_tag [run_capture_get_res_slurm_wait_finished][Error] We got something else than pid ($res[$#res]), line ".__LINE__.", run_capture_slurm called from line $line\n";
	}

	# we wait that cluster job is finnished
	while(not is_run_capture_slurm_finished($pid,$line))
	{
		sleep $SLURM_FINISHED_TIME_REQUEST;
	}
	
	# now we can get output in the dedicated cluster file

	# get the name of the process output file
	my @cluster_out_f = glob($SIGFFRID_DIR."cluster_errors/slurm*.$pid.out");

	my $nb_f = scalar(@cluster_out_f);
	if($nb_f == 1)
	{
		# we get output and store it in the array given as argument
		open(SCOF,'<',$cluster_out_f[0])or die "$prog_tag [Error] Cannot open $cluster_out_f[0] file:$!, line ".__LINE__."\n";
		while(<SCOF>)
		{
			push @$r_msg, $_;
		}
		close SCOF;
	}
	else
	{
		die "$prog_tag [run_capture_get_res_slurm_wait_finished][Error] We find $nb_f cluster output for pid:$pid, we might find only one, line ".__LINE__."\n";
	}
}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_get_res_sge_wait_finished($$$;@)

=head3 DESCRIPTION

To run linux command through capture perl lib using SGE (cluster)
when we want to get results.
Do not give back the hand until process is finished.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use sge

=item line: the line (location) where this sub is called (__LINE__)

=item r_msg: ref to an array where we store output (got from cluster out
      in fact)
      
=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 DEPENDS ON

=over

=item is_run_capture_sge_finished($$)

=item run_capture_get_res($$$;@)

=back

=cut

sub run_capture_get_res_sge_wait_finished($$$;@)
{
	my($cmd,$line,$r_msg,@error_allowed) = @_;
	run_capture_get_res_sge_prefix_wait_finished($cmd,$line,'',$r_msg,@error_allowed);
}
# ---------------------------------------------------------------------------------



# to test sge capture submission
if($b_test_run_capture_sge)
{
	# test job finished
	my $cmd = "perl sleep1.pl";
	my $pid = run_capture_sge($cmd,__LINE__);
	sleep 30;
	if(is_run_capture_sge_finished($pid,__LINE__))
	{
		print "$prog_tag $cmd finished, normal\n";
	}
	else
	{
		die "$prog_tag [Error] $cmd command with pid:$pid is not finished, it might be finished, problem, line ".__LINE__."\n";
	}

	# test job not finished
	$cmd = "perl sleep20.pl";
	print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
	$pid = run_capture_sge($cmd,__LINE__);
	print "$prog_tag ran\n";
	sleep 3;
	if(not is_run_capture_sge_finished($pid,__LINE__))
	{
		print "$prog_tag $cmd not finished, normal\n";
	}
	else
	{
		die "$prog_tag [Error] $cmd command with pid:$pid is finished, it might not be, problem, line ".__LINE__."\n";
	}
}

# to test sge capture submission
if($b_test_run_capture_get_res_sge_wait_finished)
{
	my $cmd = "echo 'Hello World'";
	my @res = ();
	run_capture_get_res_sge_wait_finished($cmd,__LINE__,\@res);
	
	print join("\n","$prog_tag [TEST][run_capture_get_res_sge_wait_finished] We must display 'Hello World':",@res,' ');
}

# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_sge_prefix_wait_finished($$$;@)

=head3 DESCRIPTION

To run linux command through capture perl lib using SGE (cluster).
Do not give back the hand until process is finished.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use sge

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script name

=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 DEPENDS ON

=over

=item is_run_capture_sge_finished($$)

=item run_capture_sge_prefix($$$;@)

=back

=cut

sub run_capture_sge_prefix_wait_finished($$$;@)
{
	my($cmd,$line,$prefix,@other_args) = @_;
	my $pid = run_capture_sge_prefix($cmd,$line,$prefix,@other_args);
	
	while(not is_run_capture_sge_finished($pid,$line)){
		sleep $SGE_FINISHED_TIME_REQUEST;
	}
}

if($b_test_run_capture_sge_prefix_wait_finished)
{
	my $test_dir = $SIGFFRID_DIR.'test_run_capture_pm/';
	my $snakefile = $test_dir.'SnakeFile_bash_echo_sge';
	my $expected_res_f = $test_dir.'snakemake_call/toto_f.txt';
	-e $expected_res_f and unlink $expected_res_f;
	my $expected_perl_res_f = $test_dir.'normal_call/toto_f_perl.txt';
	-e $expected_perl_res_f and unlink $expected_perl_res_f;
	
	-e $test_dir or die join(' ', "$prog_tag [b_test_run_capture_sge_wait_finished_force_exit_status][Error] $test_dir",
		"test_dir not found, please check, line ".__LINE__."\n");
		
	`snakemake -p -s $snakefile --cluster \"qsub -q long.q -V \" -j 2`;
	
	my $cmd = "echo toto > $expected_perl_res_f";
	print "cmd:$cmd, line ".__LINE__."\n";
	run_capture::run_capture_sge_prefix_wait_finished($cmd, __LINE__, 'echo_toto');
	print "ran\n";
}
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_slurm_prefix_wait_finished($$$;@)

=head3 DESCRIPTION

To run linux command through capture perl lib using SLURM (cluster).
Do not give back the hand until process is finished.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use slurm

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script name

=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 DEPENDS ON

=over

=item is_run_capture_slurm_finished($$)

=item run_capture_slurm_prefix($$$;@)

=back

=cut

sub run_capture_slurm_prefix_wait_finished($$$;@)
{
	my($cmd,$line,$prefix,@other_args) = @_;
	my $pid = run_capture_slurm_prefix($cmd,$line,$prefix,@other_args);
	
	while(not is_run_capture_slurm_finished($pid,$line)){
		sleep $SLURM_FINISHED_TIME_REQUEST;
	}
}

if($b_test_run_capture_slurm_prefix_wait_finished)
{
	my $test_dir = $SIGFFRID_DIR.'test_run_capture_pm/';
	my $snakefile = $test_dir.'SnakeFile_bash_echo_slurm';
	my $expected_res_f = $test_dir.'snakemake_call/toto_f.txt';
	-e $expected_res_f and unlink $expected_res_f;
	my $expected_perl_res_f = $test_dir.'normal_call/toto_f_perl.txt';
	-e $expected_perl_res_f and unlink $expected_perl_res_f;
	
	-e $test_dir or die join(' ', "$prog_tag [b_test_run_capture_slurm_wait_finished_force_exit_status][Error] $test_dir",
		"test_dir not found, please check, line ".__LINE__."\n");
		
	my $cmd = "echo toto > $expected_perl_res_f";
	print "cmd:$cmd, line ".__LINE__."\n";
	run_capture::run_capture_slurm_prefix_wait_finished($cmd, __LINE__, 'echo_toto');
	print "ran\n";
}
# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_prefix_wait_finished($$$;@)

=head3 DESCRIPTION

To run linux command through capture perl lib using SGE (cluster).
Do not give back the hand until process is finished.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use sge

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script name

=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 DEPENDS ON

=over

=item is_run_capture_finished($$)

=item run_capture_prefix($$$;@)

=back

=cut

sub run_capture_prefix_wait_finished($$$;@)
{
	my($cmd,$line,$prefix,@other_args) = @_;
	my $pid = run_capture_prefix($cmd,$line,$prefix,@other_args);
	
	while(not is_run_capture_finished($pid,$line)){
		sleep $FINISHED_TIME_REQUEST;
	}
}

if($b_test_run_capture_prefix_wait_finished)
{
	my $test_dir = $SIGFFRID_DIR.'test_run_capture_pm/';
	my $snakefile = $test_dir.'SnakeFile_bash_echo_sge';
	my $expected_res_f = $test_dir.'snakemake_call/toto_f.txt';
	-e $expected_res_f and unlink $expected_res_f;
	my $expected_perl_res_f = $test_dir.'normal_call/toto_f_perl.txt';
	-e $expected_perl_res_f and unlink $expected_perl_res_f;
	
	-e $test_dir or die join(' ', "$prog_tag [b_test_run_capture_wait_finished_force_exit_status][Error] $test_dir",
		"test_dir not found, please check, line ".__LINE__."\n");
		
	my $cmd = "echo toto > $expected_perl_res_f";
	print "$prog_tag [test_run_capture_prefix_wait_finishe] cmd:$cmd, line ".__LINE__."\n";
	run_capture::run_capture_prefix_wait_finished($cmd, __LINE__, 'echo_toto');
	print "$prog_tag [test_run_capture_prefix_wait_finishe] ran\n";
}
# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_sge_prefix_wait_finished_force_exit_status($$$;@)

=head3 DESCRIPTION

To run linux command through capture perl lib using SGE (cluster).
Do not give back the hand until process is finished.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use sge

=item line: the line (location) where this sub is called (__LINE__)

=item prefix: string inserted near the beginning of the script name

=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 DEPENDS ON

=over

=item is_run_capture_sge_finished($$)

=item run_capture_sge_prefix_force_exit_status($$$;@)

=back

=cut

sub run_capture_sge_prefix_wait_finished_force_exit_status($$$;@)
{
	my($cmd,$line,$prefix,@other_args) = @_;
	my $pid = run_capture_sge_prefix_force_exit_status($cmd,$line,$prefix,@other_args);
	
	while(not is_run_capture_sge_finished($pid,$line)){
		sleep $SGE_FINISHED_TIME_REQUEST;
	}
}

if($b_test_run_capture_sge_prefix_wait_finished_force_exit_status)
{	
	my $test_dir  = $SIGFFRID_DIR.'test_run_capture_pm/';
	my $gblock_in = $test_dir.'O157H7_and_plasmids_6_refs_merged.fasta';
	my $gblock_out = $test_dir.'O157H7_and_plasmids_6_refs_merged.fasta-gb';
	-e $gblock_out and unlink $gblock_out;

	my $cmd = "Gblocks -t=d -b5=n -b4=10 $gblock_in || true " ;
	my $prefix = "test_Gblocks";
	print "$prog_tag [test_run_capture_sge_prefix_wait_finished_force_exit_status] cmd:$cmd\n";
	run_capture_sge_prefix_wait_finished_force_exit_status($cmd,__LINE__,$prefix);
	print "$prog_tag [test_run_capture_sge_prefix_wait_finished_force_exit_status] ran\n";

}
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

=head2  NAME

run_capture_sge_wait_finished($$;@)

=head3 DESCRIPTION

To run run_capture_sge_prefix_wait_finished with an empty prefix.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use sge

=item line: the line (location) where this sub is called (__LINE__)
 
=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 DEPENDS ON

=over

=item is_run_capture_sge_finished($$)

=item run_capture_sge_prefix($$;@)

=back

=cut

sub run_capture_sge_wait_finished($$;@)
{
	my($cmd,$line,@other_args) = @_;
	my $pid = run_capture_sge_prefix($cmd,$line,'',@other_args);
	
	while(not is_run_capture_sge_finished($pid,$line)){
		sleep $SGE_FINISHED_TIME_REQUEST;
	}
}


# ---------------------------------------------------------------------------------


=head2  NAME

run_capture_force_ncbi($$;@)

=head3 DESCRIPTION

To run linux command through capture perl lib when we dont want get output,
case of connection problem with ncbi. Will retry ten times.

=head3 ARGUMENTS

=over

=item cmd:  the command line to run (string), use sge

=item line: the line (location) where this sub is called (__LINE__)
 
=item errors_allowed: list of errors allowed (and caught, integers), EXIT_ANY for all.

=back

=head3 GLOBAL VAR

number of trials:10

=cut

my $nb_trials = 10;
sub run_capture_force_ncbi($$;@)
{
	my($cmd,$line,@error_allowed) = @_;

	my $msg = '';
	my $nb_try = 0;

	while($nb_try < $nb_trials)
	{
		# run here
		try
		{
			if(@error_allowed)
			{
			  if($error_allowed[0] eq 'EXIT_ANY')
			    {
				$msg = capture(EXIT_ANY,$cmd);
			      }
			  else
			    {
			      $msg = capture([@error_allowed],$cmd);
			    }

			}
			else
			{
				$msg = capture($cmd);
			}
			$nb_try++;
		}
		catch
		{
			warn "$prog_tag [Error] cmd:$cmd failed, msg:$msg, line $line, error:$@\nWe try again\n";
			if($nb_try == $nb_trials)
			{
				die "$prog_tag [Error] cmd:$cmd failed, msg:$msg, line $line, error:$@\n10 tries failed\n";
			}
		};
		#else
		#{
		#	last;
		# }
	}
	print "$msg\n";
}
# ---------------------------------------------------------------------------------

1;
