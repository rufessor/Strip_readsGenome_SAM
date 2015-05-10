#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;     #Do not buffer output to std error

######Copyright 2015 Matthew F. Hockin Ph.D.
#All rights reserved
#The University of Utah
######mhockin@gmail.com

use File::Basename;
use Getopt::Long;

my ($I, $filePath, $prefix, $reverse, $help);
GetOptions('f=s', => \$filePath,
    'help', => \$help,       
    'rg=s', => \$I,        
    'prefix=s' => \$prefix,  
    'reverse'=> \$reverse);

if ($help == 1){
    Help();
    exit;
}

my ($SAMfile, $dir, $ext) = fileparse($filePath);
open (my $SAM_in_FH, '<', $dir.$SAMfile.$ext) or die "Cannot open $SAMfile.$ext using this path $dir\nPlease check command line input.\nLooking for GenomeMark (M,H, etc) file to edit and prefix to append to generate output\n";
open (my $newSAM_FH, '>', $dir."$prefix"."_".$SAMfile.$ext) or die "Cannot create file $SAMfile.$prefix.$ext in $dir\nPlease check command line input.\n";
Reverse($I, $prefix, $filePath) if ($reverse == 1);
Normal_mode($I, $prefix, $filePath) unless ($reverse == 1);
close $SAM_in_FH;
close $newSAM_FH;

sub Normal_mode{
    my ($I, $prefix, $filePath) = @_;
    while (my $line = <$SAM_in_FH>){
        chomp $line;
        my @lineSAM = split (/\t/, $line);
        if ($lineSAM[0] =~ m/^@/) {
            next if ($lineSAM[1] =~ m/^@SN:$I(?:MT|[XY0-9]{1,2})/);
            print $newSAM join ("\t", @lineSAM) . "\n";
            next;
        }
        next if ( $lineSAM[6] =~ /$I^(?:MT|[XY0-9]{1,2})/ );# || $lineSAM[2] =~ /^I$(?:MT|[XY0-9])/ 
        my @XA_SA = grep { $lineSAM[$_] =~ /^XA|^SA/ } (11 .. $#lineSAM);
        unless (@XA_SA){
            print $newSAM join ("\t", @lineSAM) . "\n";
            next;
        }
        my (%new_BWA, @spliceOut);
        LINE: for my $index (@XA_SA){
            my @bwa_alt_field = split (';' , $lineSAM[$index]);
            $bwa_alt_field[0] =~ s/^([XS]A:[A-Z]{1}:)//;
            my $alt_BWA_head = $1 if ($1);
            @bwa_alt_field = grep { !/^$I(?:MT|[0-9XY])/ } @bwa_alt_field;
            unless (@bwa_alt_field){
                push @spliceOut, $index;
            }else{
                $bwa_alt_field[0] = "$alt_BWA_head"."$bwa_alt_field[0]";
                my $new_field = join (';' , @bwa_alt_field).";";
                $new_BWA{$index} = $new_field;
            }
        }
        while (my ($i, $parsed_Field) = each %new_BWA){
            $lineSAM[$i] = $parsed_Field
        }
        splice (@lineSAM, $_, 1) for (@spliceOut);
        print $newSAM join ("\t" , @lineSAM) . "\n";
    }
}

sub Reverse{
    my ($I, $prefix, $filePath) = @_;
    while (my $line = <$SAM_in_FH>){
        chomp $line;
        my @lineSAM = split (/\t/, $line);
        if ($lineSAM[0] =~ m/^@/) {
            next if ($lineSAM[1] =~ m/^@SN:(?:MT|[XY0-9]{1,2})/);
            print $newSAM join ("\t", @lineSAM) . "\n";
            next;
        }
        next if ( $lineSAM[6] =~ /^(?:MT|[XY0-9]{1,2})/ );# || $lineSAM[2] =~ /^I$(?:MT|[XY0-9])/ 
        my @XA_SA = grep { $lineSAM[$_] =~ /^XA|^SA/ } (11 .. $#lineSAM);
        unless (@XA_SA){
            print $newSAM join ("\t", @lineSAM) . "\n";
            next;
        }
        my (%new_BWA, @spliceOut);
        LINE: for my $index (@XA_SA){
            my @bwa_alt_field = split (';' , $lineSAM[$index]);
            $bwa_alt_field[0] =~ s/^([XS]A:[A-Z]{1}:)//;
            my $alt_BWA_head = $1 if ($1);
            @bwa_alt_field = grep { /^$I(?:MT|[0-9XY])/ } @bwa_alt_field;
            unless (@bwa_alt_field){
                push @spliceOut, $index;
            }else{
                $bwa_alt_field[0] = "$alt_BWA_head"."$bwa_alt_field[0]";
                my $new_field = join (';' , @bwa_alt_field).";";
                $new_BWA{$index} = $new_field;
            }
        }
        while (my ($i, $parsed_Field) = each %new_BWA){
            $lineSAM[$i] = $parsed_Field
        }
        splice (@lineSAM, $_, 1) for (@spliceOut);
        print $newSAM join ("\t" , @lineSAM) . "\n";
    }

}

sub Help{
    print "This program requires output from samtools view -Lh single_genome.bed BAM_file_with_mixed_mapping.BAM > singleGenome.SAM\n";
    print "SAM/BAM files derived from BWA mapping to mixed genomic (human + mouse) index files can break things\n";
    print "This program will remove all traces of either genome and generate clean SAM files\n";
    print "This program presumes the Chr designation in the .fai used in mapping to the primary genome was [0-9XY] + MT\n";
    print "The \"Alternate Genome\" refers to the genome with a .fai where the Chr designations are as above with a prefix\n";
    print "\nThis program utilizes the --rg tag as a PREFIX in matching to 0-9XY+MT for chromosome designation\n\n";
    print "E.G. if your genome to be removed has chromosomes labeled as label1 label2... use \"label\" as --rg (read on)\n\n";
    print "Normal mode of operation removes the genome with the prefix you specify\n";
    print "Reverse mode removes the genome WITHOUT the prefix you specify\n\n";
    print "command line options use --opt=value\n\n";
    print "--f incoming .SAM file out from samtools view -Lh genome_to_keep.bed original.BAM\n\n";
    print "--rg genome flag to remove (remove genome). If human genome is 1-22,XY; Mouse is M1-M19,XY- use --rg=M to remove Mouse\n\n";
    print "--prefix Provide prefix to be appended to input file name to generate output\n\n";
    print "--reverse (NO VALUE) Use this flag to switch operation (e.g. KEEP genome tag provided in --rg=\n";
}
