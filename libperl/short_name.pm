=head1  NAME

short_name.pm

=head1 DESCRIPTION

Get file name without path. Each subprogram provide some variants.

=head1  USAGE

=over

=item short_name($file_to_get_short_name_from)

Get short name without path and WITHOUT (last) extension.

=item short_name_with_ext($file_to_get_short_name_from)

Get short name without path and WITH extension.

=item short_name_50l_max($file_to_get_short_name_from[, $max_nb_l=50])

Get short name without path and without extension, cuting
name to make it as long as max_nb_l (50 by default but
can be provided as second argument).

=item examples:

=item use lib 'libperl';

=item (...)

=item short_name_50l_max($unpaired_reads_f, 40)

=back

=cut

# **********************************************************************
# to get short name of a file, without extension
# **********************************************************************
sub short_name($)
{
	my($sh) = @_;
	# supress starting point
	$sh =~ s/^.+\///;
	# suppress compression extension
	$sh  =~ s/\.(?:gz|zip|bz)$//;
	# suppress extension	
	$sh =~ s/\.[^\.]+$//;
	return $sh;
}
# **********************************************************************
# to get short name of a file, keeping extension
# **********************************************************************
sub short_name_with_ext($){
  my($sh) = @_;
  $sh =~ s/^.+\///;
  return $sh;
}
# **********************************************************************
# to get short name of a file, without extension, limiting file name to
# 50 letters
# **********************************************************************
sub short_name_50l_max($;$)
{
	my($sh, $max_nb_l) = @_;

	defined $max_nb_l or $max_nb_l = 50;
	
	# supress starting point
	$sh =~ s/^.+\///;
	# suppress extension
	$sh =~ s/\.[^\.]+$//;

	if(length($sh) > $max_nb_l)
	{
		$sh = substr($sh, 0, $max_nb_l);
	}
	return $sh;
}
# **********************************************************************
1;
