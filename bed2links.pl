#!/usr/bin/perl
use strict;
use warnings;
use autodie;

my @links;
open my $BED, "<", $ARGV[0];
while (<$BED>) {
    chomp;
    my @tmp = split( "\t", $_ );
    if ( $tmp[5] eq "F" ) {
        if ( $tmp[4] != 0 ) {
            ${ $links[ $tmp[3] - 1 ] }[ $tmp[4] - 1 ] =
"Chr$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]-$tmp[4]\[$tmp[0](+):$tmp[1]-$tmp[2]\]$tmp[6]\t1\t+";
        }
    }
    else {
        if ( $tmp[4] != 0 ) {
            if ( $tmp[4] == 1 ) {
                ${ $links[ $tmp[3] - 1 ] }[ $tmp[4] - 1 ] =
"Chr$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]-$tmp[4]\[$tmp[0](+):$tmp[1]-$tmp[2]\]$tmp[6]\t1\t+";
            }
            else {
                ${ $links[ $tmp[3] - 1 ] }[ $tmp[4] - 1 ] =
"Chr$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]-$tmp[4]\[$tmp[0](-):$tmp[1]-$tmp[2]\]$tmp[6]\t1\t-";
            }
        }
    }
}
close($BED);

foreach my $key ( 0 .. $#links ) {
    if ( defined( $links[$key] ) ) {
        print( join( "\n", @{ $links[$key] } ) );
        print("\n");
    }
}
__END__
