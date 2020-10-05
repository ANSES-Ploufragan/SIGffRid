#!/usr/bin/perl
$|++;

=head1  NAME

quicksortnummultiarray.pm

=head1 DESCRIPTION

Subprograms dedicated to sort NUMERICALLY on range defined by the two first arguments
the first array. Sort all other arrays according to this first array.
Composed of 6 subprograms.

=head1  USAGE

=over

=item [partition]     			makes two indexes browsing from down and up of range until upper index
								has value lower than lower index: exchange values then.
								Increasing order.

=item [partition_reverse]		makes two indexes browsing from down and up of range until upper index
								has value lower than lower index: exchange values then.
								Decreasing order.
								
=item [quicksort_recurse] 		Deal recursion of the 'quick sort' algorithm
								for increasing order
								
=item [quicksort_recurse_reverse] Deal recursion of the 'quick sort' algorithm
                                for decreasing order

=item quicksortnum      	The main program to be used externally to sort
							in increasing order

=item quicksortnum_reverse 	The main program to be used externally to sort
							in decreasing order
                                            
=back

=cut

# **********************************************************************
=head2  NAME

partition

=head2 DESCRIPTION

Make two indexes browsing from down and up of range until upper index
has value lower than lower index: exchange values then.
Increasing order.

=head2  USAGE

=over

=item partition(	$ref_to_array_sorted,
					$first_index_to_sort,
					$last_index_to_sort,
					@list_of_refs_to_other_arrays_to_sort)

=back

=cut

sub partition {
    my ( $array, $first, $last , @arrays) = @_;

    my $i = $first;
    my $j = $last - 1;
    my $pivot = $array->[ $last ];

 SCAN: {
        do {
            # $first <= $i <= $j <= $last - 1
            # Point 1.

            # Move $i as far as possible.
            while ( $array->[ $i ] <= $pivot ) {  
                $i++;
		if(not defined $pivot){
		  die "pivot not defined (index $last in the array)\n";
		}
		if(not defined $array->[ $i ]){
		  die "array of $i  not defined\n";
		}
                last SCAN if $j < $i;
            }

            # Move $j as far as possible.
            while ( $array->[ $j ] >= $pivot ) {
                $j--;
                last SCAN if $j < $i;
            }

# $i and $j did not cross over, so swap a low and a high value.
            @$array[ $j, $i ] = @$array[ $i, $j ];
	    foreach(@arrays){
	      @$_[ $j, $i ] = @$_[ $i, $j ];
	    }

        } while ( --$j >= ++$i );
    }
    # $first - 1 <= $j < $i <= $last
    # Point 2.

# Swap the pivot with the first larger
# element (if there is one).
    if ( $i < $last ) {
        @$array[ $last, $i ] = @$array[ $i, $last ];
	foreach(@arrays){
	  @$_[ $last, $i ] = @$_[ $i, $last ];
	}
        ++$i;
    }


    ++$i while $i <= $last  && $array->[ $i ] == $pivot;
    --$j while $j >= $first && $array->[ $j ] == $pivot;

    # Point 3.

    return ( $i, $j );   # The new bounds exclude the middle.
}
# **********************************************************************

# **********************************************************************
=head2  NAME

partition_reverse

=head2 DESCRIPTION

Make two indexes browsing from down and up of range until upper index
has value lower than lower index: exchange values then.
Decreasing order.

=head2  USAGE

=over

=item partition_reverse(	$ref_to_array_sorted,
							$first_index_to_sort,
							$last_index_to_sort,
							@list_of_refs_to_other_arrays_to_sort)

=back

=cut

# for decreasing order
sub partition_reverse {
    my ( $array, $first, $last , @arrays) = @_;

    my $i = $first;
    my $j = $last - 1;
    my $pivot = $array->[ $last ];

 SCAN: {
        do {
            # $first <= $i <= $j <= $last - 1
            # Point 1.

            # Move $i as far as possible.
            while ( $array->[ $i ] >= $pivot ) {  
                $i++;
		if(not defined $pivot){
		  die "pivot not defined (index $last in the array)\n";
		}
		if(not defined $array->[ $i ]){
		  die "array of $i  not defined\n";
		}
                last SCAN if $j < $i;
            }

            # Move $j as far as possible.
            while ( $array->[ $j ] <= $pivot ) {
                $j--;
                last SCAN if $j < $i;
            }

# $i and $j did not cross over, so swap a low and a high value.
            @$array[ $j, $i ] = @$array[ $i, $j ];
	    foreach(@arrays){
	      @$_[ $j, $i ] = @$_[ $i, $j ];
	    }

        } while ( --$j >= ++$i );
    }
    # $first - 1 <= $j < $i <= $last
    # Point 2.

# Swap the pivot with the first larger
# element (if there is one).
    if ( $i < $last ) {
        @$array[ $last, $i ] = @$array[ $i, $last ];
	foreach(@arrays){
	  @$_[ $last, $i ] = @$_[ $i, $last ];
	}
        ++$i;
    }


    ++$i while $i <= $last  && $array->[ $i ] == $pivot;
    --$j while $j >= $first && $array->[ $j ] == $pivot;

    # Point 3.

    return ( $i, $j );   # The new bounds exclude the middle.
}
# **********************************************************************

# **********************************************************************
=head2  NAME

quicksort_recurse

=head2 DESCRIPTION

Deal recursion of the 'quick sort' algorithm for increasing order.

=head2  USAGE

=over

=item quicksort_recurse(	$ref_to_array_sorted,
							$first_index_to_sort,
							$last_index_to_sort,
							@list_of_refs_to_other_arrays_to_sort)

=back

=cut

sub quicksort_recurse {
    my ( $array, $first, $last , @arrays) = @_;

    if ( $last > $first ) {
        my ( $first_of_last, $last_of_first ) =
                                partition( $array, $first, $last, @arrays );

        local $^W = 0; # Silence deep recursion warning.
        quicksort_recurse($array, $first,         $last_of_first, @arrays);
        quicksort_recurse($array, $first_of_last, $last, @arrays);
    }
}
# **********************************************************************

# **********************************************************************
=head2  NAME

quicksort_recurse_reverse

=head2 DESCRIPTION

Deal recursion of the 'quick sort' algorithm for decreasing order.

=head2  USAGE

=over

=item quicksort_recurse_reverse(	$ref_to_array_sorted,
									$first_index_to_sort,
									$last_index_to_sort,
									@list_of_refs_to_other_arrays_to_sort)

=back

=cut

sub quicksort_recurse_reverse {
    my ( $array, $first, $last , @arrays) = @_;

    if ( $last > $first ) {
        my ( $first_of_last, $last_of_first ) =
                                partition_reverse( $array, $first, $last, @arrays );

        local $^W = 0; # Silence deep recursion warning.
        quicksort_recurse_reverse($array, $first,         $last_of_first, @arrays);
        quicksort_recurse_reverse($array, $first_of_last, $last, @arrays);
    }
}
# **********************************************************************

# **********************************************************************
# increasing order
=head2  NAME

quicksortnum

=head2 DESCRIPTION

The main program to be used externally to sort in increasing order.

=head2  USAGE

=over

=item quicksortnum(	
									$first_index_to_sort,
									$last_index_to_sort,
									$ref_to_array_sorted,
									@list_of_refs_to_other_arrays_to_sort)

=back

=cut

sub quicksortnum ($$$@){
# ibeg, iend, reftoarray, otherarraystosort
# The recursive version is bad with BIG lists
# because the function call stack gets REALLY deep.
    quicksort_recurse($_[ 2 ], $_[ 0 ], $_[ 1 ], @_[3..$#_]);
}
# **********************************************************************

# **********************************************************************
# decreasing order
=head2  NAME

quicksortnum_reverse

=head2 DESCRIPTION

The main program to be used externally to sort in decreasing order.

=head2  USAGE

=over

=item quicksortnum(	
									$first_index_to_sort,
									$last_index_to_sort,
									$ref_to_array_sorted,
									@list_of_refs_to_other_arrays_to_sort)

=back

=cut
sub quicksortnum_reverse ($$$@){
# ibeg, iend, reftoarray, otherarraystosort
# The recursive version is bad with BIG lists
# because the function call stack gets REALLY deep.
    quicksort_recurse_reverse($_[ 2 ], $_[ 0 ], $_[ 1 ], @_[3..$#_]);
}
# **********************************************************************

1;

# If you expect that many of your keys will be the same,
# try adding this before the <LITERAL>return</LITERAL> in
# <LITERAL>partition()</LITERAL>:
#
# Extend the middle partition as much as possible.
#
# ++$i while $i <= $last  && $array->[ $i ] eq $pivot;
# --$j while $j >= $first && $array->[ $j ] eq $pivot;



# TEST ************************************************
# @array = qw(40 2 68 6 20 50 14 12 12 12 49);
# @array2 = qw(1 2 3 4 5 6 7 8 9 10 11);

# print "BEFORE:\n";
# print "@array\n";
# print "@array2\n";

# quicksort( \@array, 0, $#array, \@array2 );

# print "AFTER:\n";
# print "@array\n";
# print "@array2\n";
