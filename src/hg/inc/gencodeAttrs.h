/* gencodeAttrs.h was originally generated by the autoSql program, which also 
 * generated gencodeAttrs.c and gencodeAttrs.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef GENCODEATTRS_H
#define GENCODEATTRS_H

#define GENCODEATTRS_NUM_COLS 14

extern char *gencodeAttrsCommaSepFieldNames;

struct gencodeAttrs
/* Basic set of attributes associated with all Gencode transcripts. */
    {
    struct gencodeAttrs *next;  /* Next in singly linked list. */
    char *geneId;	/* Gene identifier */
    char *geneName;	/* Gene name */
    char *geneType;	/* BioType of gene */
    char *unused1;	/* unused (was geneStatus in wgGencode tracks) */
    char *transcriptId;	/* Transcript identifier */
    char *transcriptName;	/* Transcript name */
    char *transcriptType;	/* BioType of transcript */
    char *unused2;	/* unused (was transcriptStatus in wgGencode tracks) */
    char *unused3;	/* unused (was havanaGeneId in wgGencode tracks) */
    char *unused4;	/* unused (was havanaTranscriptId in wgGencode tracks) */
    char *ccdsId;	/* CCDS identifier if transcript is in CCDS */
    int level;	/* GENCODE level: 1 = experimental confirmed, 2 = manual, 3 = automated */
    char *transcriptClass;	/* high level type of transcript */
    char *proteinId;	/* Protein identifier (not loaded on many older versions of GENCODE) */
    };

void gencodeAttrsStaticLoad(char **row, struct gencodeAttrs *ret);
/* Load a row from gencodeAttrs table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct gencodeAttrs *gencodeAttrsLoad(char **row);
/* Load a gencodeAttrs from row fetched with select * from gencodeAttrs
 * from database.  Dispose of this with gencodeAttrsFree(). */

struct gencodeAttrs *gencodeAttrsLoadAll(char *fileName);
/* Load all gencodeAttrs from whitespace-separated file.
 * Dispose of this with gencodeAttrsFreeList(). */

struct gencodeAttrs *gencodeAttrsLoadAllByChar(char *fileName, char chopper);
/* Load all gencodeAttrs from chopper separated file.
 * Dispose of this with gencodeAttrsFreeList(). */

#define gencodeAttrsLoadAllByTab(a) gencodeAttrsLoadAllByChar(a, '\t');
/* Load all gencodeAttrs from tab separated file.
 * Dispose of this with gencodeAttrsFreeList(). */

struct gencodeAttrs *gencodeAttrsCommaIn(char **pS, struct gencodeAttrs *ret);
/* Create a gencodeAttrs out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new gencodeAttrs */

void gencodeAttrsFree(struct gencodeAttrs **pEl);
/* Free a single dynamically allocated gencodeAttrs such as created
 * with gencodeAttrsLoad(). */

void gencodeAttrsFreeList(struct gencodeAttrs **pList);
/* Free a list of dynamically allocated gencodeAttrs's */

void gencodeAttrsOutput(struct gencodeAttrs *el, FILE *f, char sep, char lastSep);
/* Print out gencodeAttrs.  Separate fields with sep. Follow last field with lastSep. */

#define gencodeAttrsTabOut(el,f) gencodeAttrsOutput(el,f,'\t','\n');
/* Print out gencodeAttrs as a line in a tab-separated file. */

#define gencodeAttrsCommaOut(el,f) gencodeAttrsOutput(el,f,',',',');
/* Print out gencodeAttrs as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* GENCODEATTRS_H */

