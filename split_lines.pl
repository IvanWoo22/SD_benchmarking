#!/usr/bin/env perl
use strict;
use warnings;

open my $input_fh, '<', $ARGV[0] or die "$!";
my @lines = <$input_fh>;
close $input_fh;

sub replace_chromosome_names {
    my ($line) = @_;
    $line =~ s/Chr//g;
    return $line;
}

open my $output_fh, '>', $ARGV[1] or die "$!";
foreach my $line (@lines) {
    chomp $line;
    my @fields = split /\t/, $line;
    for my $i ( 0 .. $#fields - 1 ) {
        for my $j ( $i + 1 .. $#fields ) {
            my ( $chr1, $chr2 ) =
              map { s/\([+-]\):(\d+)-(\d+)/ $1 $2/g; $_ } @fields[ $i, $j ];
            $chr1 = replace_chromosome_names($chr1);
            $chr2 = replace_chromosome_names($chr2);
            print $output_fh "$chr1 $chr2\n";
        }
    }
}
close $output_fh;

__END__
