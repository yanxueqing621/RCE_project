#!/usr/bin/perl
use IO::All;
use Modern::Perl;
use Data::Dumper;

my $anno_file = 'gencode.v44.annotation_longest_cds.gtf.anno';
my $io_anno = io("$anno_file")->chomp;
$io_anno->getline;

my %anno;
while (my $line = $io_anno->getline) {
    my @cols = split /\t/, $line;
    #    say Dumper \@cols;
    $anno{$cols[0]}{l_cds}  = $cols[3];
    $anno{$cols[0]}{l_utr3}  = $cols[4];
    $anno{$cols[0]}{s_utr5}  = 1;
    $anno{$cols[0]}{e_utr5}  = $cols[2];
    $anno{$cols[0]}{s_cds}  = $cols[2]+1;
    $anno{$cols[0]}{e_cds}  = $cols[2]+$cols[3]-3;
    $anno{$cols[0]}{e_stop}  = $cols[2]+$cols[3];
    $anno{$cols[0]}{s_utr3}  = $cols[2] + $cols[3] + 1;
    $anno{$cols[0]}{e_utr3}  = $cols[1];
    $anno{$cols[0]}{stop} = $cols[5];
    $anno{$cols[0]}{symbol} = $cols[-1];
}

#summary_sample("MTGA_psite_29nt_test.xls");
summary_sample(shift);

sub summary_sample {
    my $psite_file = shift;
    my $io_psite = io("$psite_file")->chomp;
    #my $io_out2 = io("$psite_file.region");
    $io_psite->getline;
    my %genes;
    my %samples;
    while (defined (my $line = $io_psite->getline)) {
        #chop $line;  #remove
        my @cols = split /\s+/, $line;
        next if not defined $anno{$cols[0]};
        my $stop_codon = $anno{$cols[0]}{stop};
        if ($cols[1] <= $anno{$cols[0]}{e_utr5}) {
            $samples{$psite_file}{$stop_codon}{utr5}++;
            $genes{$cols[0]}{utr5}++;
            #$io_out2->println("$line\t5utr");
        } elsif ($cols[1] <=$anno{$cols[0]}{e_cds}) {
            $samples{$psite_file}{$stop_codon}{cds}++;
            $genes{$cols[0]}{cds}++;
            #$io_out2->println("$line\tcds");
        } elsif ($cols[1] <= $anno{$cols[0]}{e_stop}) {
            next;  #NTC  stop codon skip
        } elsif ($cols[1] <= $anno{$cols[0]}{e_utr3} and $cols[1] >= $anno{$cols[0]}{e_cds} + 6) {
            $samples{$psite_file}{$stop_codon}{utr3}++;
            $genes{$cols[0]}{utr3}++;
            #$io_out2->println("$line\t3utr");
        } else {
            next;
        }
        $samples{$psite_file}{$stop_codon}{trans}++;
        $genes{$cols[0]}{trans}++;
    }
    for my $stop_codon (qw/TGA TAA TAG/) {
        $samples{$psite_file}{$stop_codon}{ratio} = $samples{$psite_file}{$stop_codon}{utr3} / $samples{$psite_file}{$stop_codon}{cds};
    }
    say Dumper \%samples;

    # calculate genes RRTS score
    my $io_gene_rrts = io("$psite_file.gene.RRTS");
    $io_gene_rrts->println("gid\tgene_symbol\tgene_stop_codon\tcds_len\tutr3_len\tcds_count\tutr3_count\trrts_score");
    for my $gid (keys %genes) {
       next if not defined $genes{$gid}{cds};
       my $symbol = $anno{$gid}{symbol};
       my $stop_codon = $anno{$gid}{stop};
       my $l_cds =  $anno{$gid}{l_cds};
       my $l_utr3 =  $anno{$gid}{l_utr3};
       my $c_cds = $genes{$gid}{cds};
       my $c_utr3 = defined $genes{$gid}{utr3} ? $genes{$gid}{utr3} : 0;
       my $line = join "\t", $gid, $symbol, $stop_codon, $l_cds, $l_utr3, $c_cds, $c_utr3, ($c_utr3/$l_utr3)/($c_cds/$l_cds);
       $io_gene_rrts->println($line); #if $stop_codon eq 'TGA';
    }
}

