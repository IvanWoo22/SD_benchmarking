#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

my @chr_set;
foreach ( 1 .. 5 ) {
    $chr_set[$_] = AlignDB::IntSpan->new;
}

open GFF, "<", $ARGV[0];
while (<GFF>) {
    chomp;
    my @tmp   = split( "\t", $_ );
    my $chr   = $tmp[0];
    my $start = $tmp[1];
    my $end   = $tmp[2];
    $chr_set[$chr]->AlignDB::IntSpan::add_pair( $start, $end );
}
close GFF;

open SEG, "<", $ARGV[1];
while (<SEG>) {
    chomp;
    my @tmp = split( "\t", $_ );
    my $set = AlignDB::IntSpan->new;
    $set->add_pair( $tmp[1], $tmp[2] );
    my $intsec = $set->intersect( $chr_set[ $tmp[0] ] );

    if ( $intsec->cardinality >= 5 ) {
        my $rl = $intsec->as_string;
        print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$rl\n";
    }
}
close SEG;

__END__