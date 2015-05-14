# Strip_readsGenome_SAM
HTS reads dervied from mixed genome samples must be mapped to combined genomic index files (using for instance BWA)
Although this is required for correct read placement it can be problematic for downstream processing

This script is a catch all to work through some of the SAM file manipulations requied to deal with these issues.

In general- this script can accomplish the following things while either writing a new .SAM file prepended with the --prefix=NAME or through streaming output to the command line by invoking --stream.

It has the following dependicies-
File::Basename
Getopt::Long

the command line options are 
--prefix=<> (provide name to be prepended to input file to generate unique output file)
--rg=<>
--help (display abbreviated help)
--reverse (switch mode from removal of rg specified tags to removal of tags NOT corresponding to rg)
--insert (switch mode entirely, now we prepend rg to the reference seq name- refSeq name must correspond to [0-9]|X|Y|MT)

This script is built to edit BWA output and examines the alt SAM fields BWA generates.  This script operates on SAM files from paired end or single ended reads- it operates on the RNAME and RNEXT fields. In addition it parses the BWA generated alt fields XA|SA in the same manner as it does for the RNAME and RNEXT fields.  

For example, if you have mapped to a combined Human and Mouse index where in the Mouse SN fields have been altered from the typical [0-9]|X\Y|MT to instead include a prefix "M"- e.g. M1, MMT, MX, MY for the sequence names, and the human are represented normally, e.g. 1,2,X,Y,MT- and you wish to retain ONLY reads that map to the human genome.

perl Remove_genomeV3.pl --rg=M --prefix=StripGRCm /path/to/sam/MYfile.SAM 

Will generate a new file- StripGRCm_MYfile.SAM in which any read or read with a mate mapping to any reference sequence name with a M as a prefix has been removed.  In addition, any alt mapping fields with reference to the same will be removed.  Any alt mapping fields with reference to the human genome will be retained.  In the instance that the only XA|SA mapping is to Mouse sequence names, the entire XA|SA field is removed.

AS SUCH- reads that themselves map to the human genome, but for which their mate mapped to mouse, will be removed entirely- this behavior is intended.

If instead you wanted to retain the mouse mapping and remove the human- just invoke the --reverse flag (which does not take an argument)

Finally, if you wish to modify a SAM file to INCLUDE a prefix to the standard [0-9]|X|Y|MT designation, just use --insert and the rg tag will be prepended to every instance of the reference sequence name.  This option does NOT work on mixed files- you must have a clean file with mapping to reference sequence names described as [0-9|X|Y|MT.

