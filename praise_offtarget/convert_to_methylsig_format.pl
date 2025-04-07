#!/usr/bin/perl
use IO::All;
use Modern::Perl;
use Data::Dumper;

my ($infile, $outfile) = @ARGV;
my $io_in = io("$infile")->chomp;
my $io_out = io("$outfile");

$io_in->getline;
while (defined (my $line = $io_in->getline)) {
    my ($gid, $pos, $total, $deletion,$ratio ) = split /,/, $line;
        my ($pos1, $pos2);
        if ($pos =~/-/) {
                ($pos1, $pos2) = $pos =~/(\d+)-(\d+)/;
                $pos2 = $pos2 + 1;
        } else {
                $pos1 = $pos;
                $pos2 = $pos + 1;
        }

        my $undeletion = $total - $deletion;
        my $newline = join "\t", $gid, $pos1, $pos2, $ratio, $deletion, $undeletion;
        $io_out->println($newline);
}

