#!/usr/bin/env perl
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $VERSION = "0.0.1";
my $quiet = 0;
my $EXE = basename($0);
my $seqtype = "paired";
my ($read1, $read2, $isolate);
&GetOptions(
    "t=s" =>\$seqtype,
    "1=s" =>\$read1,
    "2=s" =>\$read2,
    "s=s" =>\$isolate
    );

($seqtype && $read1) || 
    die "Name:\n".
    "by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  generated stats for the input seq files\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -t <paired|single>\n".
    "  -1 <read1 file> \n".
    "  -2 <read2 file> \n" .
    "  -s <isolate name> \n";


my @seq_files = ();
push(@seq_files, $read1);
push(@seq_files, $read2) if $seqtype eq "paired";


    
    #min_len: 35; max_len: 151; avg_len: 147.63; 6 distinct quality values
    #POS     #bases  %A      %C      %G      %T      %N      avgQ    errQ    %low    %high
    #ALL     134985892       21.7    28.2    28.6    21.5    0.0     32.1    22.9    11.3 >
    #1       914379  19.4    23.9    43.8    13.0    0.0     30.9    25.7    4.3     95.7
    #2       914379  18.3    28.6    18.0    35.1    0.0     31.1    26.3    3.6     96.4
    #....
    #150     785526  22.8    26.1    31.3    19.8    0.0     27.3    19.8    23.8    76.2
    #151     530460  31.4    0.0     44.7    23.9    0.0     23.9    17.9    38.3    61.7

    my $filename_str = join(" ", @seq_files);
    msg("processed data from @seq_files");
    
    my %stat = ();
    my $cmd = "cat $filename_str  | seqtk fqchk -q0 -";
    msg("running command: $cmd");
    open my $IN, '-|', $cmd or err("could not run command: $cmd");

    while (<$IN>) {
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
    msg("processed", $stat{num_reads}, "reads from $read1 and $read2 dataset.");
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
    

  


sub msg { print STDERR "[$EXE] @_\n" unless $quiet; }
sub err { $quiet=0; msg("ERROR:", @_); exit(-1); }

