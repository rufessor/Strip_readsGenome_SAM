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

my ($I, $prefix);
my $reverse = 0;
my $help = 0;
my $insert = 0 ;
my $stream = 0;
GetOptions('help', => \$help,       
    'rg=s', => \$I,        
    'prefix=s' => \$prefix,  
    'reverse' => \$reverse,
    'insert' => \$insert,
    'stream' => \$stream)
    or die("Error in command line arguments\n");

if ($help == 1){
    Help();
    exit;
}
die "You MUST include --prefix=NAME to output to a new file- if you want to stream output use --stream without invoking --prefix= \n"; if ($stream == 0 && !$prefix)
my $filePath = shift;
my ($SAMfile, $dir, $ext) = fileparse($filePath);
open (my $SAM_in_FH, '<', $dir.$SAMfile.$ext) or die "Cannot open your input SAM file $SAMfile.$ext using this path $dir\n\n";
open (my $newSAM_FH, '>', $dir."$prefix"."_".$SAMfile.$ext) or die "Cannot create file $prefix $SAMfile.$ext in $dir\n\n" if ($stream == 0);

Insert_mode($I, $prefix, $filePath, $newSAM_FH) if ($insert == 1 && $stream == 0);
Insert_mode($I, $prefix, $filePath) if ($insert == 1 && $stream == 1);
Reverse($I, $prefix, $filePath, $newSAM_FH) if ($reverse == 1 && $stream == 0);
Reverse($I, $prefix, $filePath) if ($reverse == 1 && $stream == 1);
Normal_mode($I, $prefix, $filePath, $newSAM_FH) if ($reverse == 0 && $stream == 0);
Normal_mode($I, $prefix, $filePath) if ($reverse == 0 && $stream == 1);
close $SAM_in_FH;
close $newSAM_FH (if $stream == 0);

sub Normal_mode{
    my ($I, $prefix, $filePath, $newSAM) = @_;
    while (my $line = <$SAM_in_FH>){
        chomp $line;
        my @lineSAM = split (/\t/, $line);
        if ($lineSAM[0] =~ m/^@/) {
            next if ($lineSAM[1] =~ m/^SN:$I(?:MT|[XY0-9]{1,2})/);
            print $newSAM join ("\t", @lineSAM) . "\n";
            next;
        }
        next if ( $lineSAM[6] =~ /^$I(?:MT|[XY0-9]{1,2})/ );# || $lineSAM[2] =~ /^I$(?:MT|[XY0-9])/ 
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
    my ($I, $prefix, $filePath, $newSAM) = @_;
    while (my $line = <$SAM_in_FH>){
        chomp $line;
        my @lineSAM = split (/\t/, $line);
        if ($lineSAM[0] =~ m/^@/) {
            next if ($lineSAM[1] =~ m/^SN:(?:MT|[XY0-9]{1,2})/);
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

#XA:Z:19,-7734532,91M1I33M,8;4,-74372427,91M1I33M,8;13,+46322743,33M1I91M,8;4,-71972611,91M1I33M,9;9,+84442877,33M1I91M,9;
sub Insert_mode{
    print STDERR "INSERT MODE \n";
    my ($I, $prefix, $filePath, $newSAM) = @_;
    while (my $line = <$SAM_in_FH>){
        chomp $line;
        my @lineSAM = split (/\t/, $line);
        if ($lineSAM[0] =~ m/^@/) {
            if ($lineSAM[1] =~ /^SN:/){
                $lineSAM[1] =~ /^SN:(MT|[XY0-9]{1,2})/;
                $lineSAM[1] = "SN:".$I.$1;
                $stream == 1 ? print join ("\t" , @lineSAM) . "\n" : print $newSAM join ("\t" , @lineSAM) . "\n";
                next;
            } 
            $stream == 1 ? print join ("\t" , @lineSAM) . "\n" : print $newSAM join ("\t" , @lineSAM) . "\n"; 
            next;
        }
        $lineSAM[2] = $I.$lineSAM[2];
        if ($lineSAM[6] ne "="){
            $lineSAM[6] = $I.$lineSAM[6];
        }
        my @XA_SA = grep { $lineSAM[$_] =~ /^[XS]A/ } (11 .. $#lineSAM);
        unless (@XA_SA){
            $stream == 1 ? print join ("\t" , @lineSAM) . "\n" : print $newSAM join ("\t" , @lineSAM) . "\n";
            next;
        }
        my (%new_BWA, @spliceOut);
        LINE: for my $index (@XA_SA){
            my @bwa_alt_field = split (';' , $lineSAM[$index]);
            $bwa_alt_field[0] =~ s/^([XS]A:[A-Z]{1}:)//;
            my $alt_BWA_head = $1;
            my @new_altField;
            push @new_altField, $I.$_ for (@bwa_alt_field);   
            #@bwa_alt_field = grep { !/^$I(?:MT|[0-9XY])/ } @bwa_alt_field;
            unless (@bwa_alt_field){
                push @spliceOut, $index;
            }else{
                $new_altField[0] = "$alt_BWA_head"."$new_altField[0]";
                my $new_field = join (';' , @new_altField).";";
                $new_BWA{$index} = $new_field;
            }
        }
        while (my ($i, $parsed_Field) = each %new_BWA){
            $lineSAM[$i] = $parsed_Field;
        }
        splice (@lineSAM, $_, 1) for (@spliceOut);
        $stream == 1 ? print join ("\t" , @lineSAM) . "\n" : print $newSAM join ("\t" , @lineSAM) . "\n";
    }
}

sub Help{
    print "This program requires output from samtools view -Lh single_genome.bed BAM_file_with_mixed_mapping.BAM > singleGenome.SAM\n";
    print "SAM/BAM files derived from BWA mapping to mixed genomic (human + mouse) index files can break things\n";
    print "This program will remove all traces of either genome and generate clean SAM files\n";
    print "Normal mode operation presumes the Chr designation in the .fai used in mapping to the primary genome was [0-9XY] + MT\n";
    print "The \"Alternate Genome\" refers to the genome with a .fai where the Chr designations are as above with a prefix\n";
    print "\nThis program utilizes the --rg tag as a PREFIX in matching to 0-9XY+MT for chromosome designation\n\n";
    print "E.G. if your genome to be removed has chromosomes labeled as label1 label2... use \"label\" as --rg (read on)\n\n";
    print "Normal mode of operation removes the genome with the prefix you specify\n";
    print "Reverse mode removes the genome WITHOUT the prefix you specify\n\n";
    print "Insert Mode adds a prefix to all fields with a genome designation in the input SAM field.  It only works if your SAM file uses [0-9XY],MT\n\n";
    print "command line options use --opt=value\n\n";
    print "specify incoming .SAM file on command line without option as relative or full path\n\n";
    print "--rg genome flag to remove (remove genome). If human genome is 1-22,XY; Mouse is M1-M19,XY- use --rg=M to remove Mouse\n\n";
    print "--prefix Provide prefix to be appended to input file name to generate output\n\n";
    print "--reverse (NO VALUE) Use this flag to switch operation (e.g. KEEP only SAM lines with primairy or mate pair mapping to --rg also discards XA/SA alt fields without genome tag provided in --rg=\n";
    print "--insert (NO VALUE) Use this flag to switch operation (e.g. INSERT the --rg= tag as a prefix to all chr designations.\n\n";
    print "--stream (NO VALUE) coverts operation from file writing (output --prefix=NAME yeilds NAME_INfile as new file) to streaming\n\n";
    print "if using --stream any input to --prefix will be disregarded and NO output file written\n\n";
}
