#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';

my $seqtype = "paired";
my ($input, $isolate);
&GetOptions(
    "t=s" =>\$seqtype,
    "i=s" =>\$input,
    "s=s" =>\$isolate
    );

($seqtype && $input) || 
    die "Name:\n".
    "by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  generated stats for the input seq files\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -t <paired|single>\n".
    "  -i <output from seqkit fqchk> \n".
    "  -s <isolate name> \n";

#min_len: 35; max_len: 151; avg_len: 147.63; 6 distinct quality values
#POS     #bases  %A      %C      %G      %T      %N      avgQ    errQ    %low    %high
#ALL     134985892       21.7    28.2    28.6    21.5    0.0     32.1    22.9    11.3 >
#1       914379  19.4    23.9    43.8    13.0    0.0     30.9    25.7    4.3     95.7
#2       914379  18.3    28.6    18.0    35.1    0.0     31.1    26.3    3.6     96.4
#....
#150     785526  22.8    26.1    31.3    19.8    0.0     27.3    19.8    23.8    76.2
#151     530460  31.4    0.0     44.7    23.9    0.0     23.9    17.9    38.3    61.7

my %stat = ();
open(IN, '<', $input) or die "Could not opne $input to read, $!\n";
while (<IN>) {
    if (m/^min_len/) {
        s/\s//g;
        for my $pair (split m';') {
        my($k,$v) = split m':', $pair;
        $stat{$k} = $v if $v;
        }
    }
    elsif (m/^ALL/) {
        my @x = split ' ';
        $stat{total_bp} = $x[1];
        $stat{gee_cee} = $x[3] + $x[4];
        $stat{avg_qual} = $x[7];
        $stat{err_qual} = $x[8];
        $stat{ambig_bp_pc} = $x[6];
    }
    elsif (m/^1\s+(\d+)\b/) {
        $stat{num_reads} = $1;
    }
}
my @out = ();
    
push(@out, $isolate);
if($seqtype eq "paired"){
    push(@out, int($stat{num_reads}/2));
}
else{
    push(@out, $stat{num_reads});
}
push(@out,  $stat{total_bp});
push(@out,  $stat{gee_cee});
push(@out,  $stat{min_len});
push(@out,  int($stat{avg_len}));
push(@out,  $stat{max_len});
push(@out, $stat{avg_qual});
push(@out, $stat{err_qual});
push(@out, $stat{ambig_bp_pc});
#print join(",", @out), "\n";

my @header = ("Isolate", "Reads", "Yield", "GeeCee", "MinLen", "AvgLen", "MaxLen", "AvgQual", "ErrQual", "Ambiguous");

print join("\t", @header), "\n";
print join("\t", @out), "\n";
