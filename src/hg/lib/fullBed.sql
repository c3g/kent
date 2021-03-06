# fullBed.sql was originally generated by the autoSql program, which also 
# generated fullBed.c and fullBed.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Browser extensible data
CREATE TABLE bed (
    bin smallint not null,              # Bin number for browser speedup
    chrom varchar(255) not null,	# chromosome or FPC contig
    chromStart int unsigned not null,	# Start position in chromosome
    chromEnd int unsigned not null,	# End position in chromosome
    name varchar(255) not null,	        # Name of item
    score uint not null,	        # Score from 0-1000
    strand char(1) not null,	        # + or -
    thickStart int unsigned not null,	# Start of where display should be thick (start codon)
    thickEnd int unsigned not null,	# End of where display should be thick (stop codon)
    reserved int unsigned not null,	# Used as itemRgb as of 2004-11-22
    blockCount int not null,	        # Number of blocks
    blockSizes longblob not null,	# Comma separated list of block sizes
    chromStarts longblob not null,	# Start positions relative to chromStart
              #Indices
    INDEX(chrom(8),bin),
    INDEX(name(20))
);
