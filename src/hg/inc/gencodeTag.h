/* gencodeTag.h was originally generated by the autoSql program, which also 
 * generated gencodeTag.c and gencodeTag.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef GENCODETAG_H
#define GENCODETAG_H

#define GENCODETAG_NUM_COLS 2

extern char *gencodeTagCommaSepFieldNames;

struct gencodeTag
/* Tags associated with GENCODE transcripts. */
    {
    struct gencodeTag *next;  /* Next in singly linked list. */
    char *transcriptId;	/* GENCODE transcript identifier */
    char *tag;	/* symbolic tag */
    };

void gencodeTagStaticLoad(char **row, struct gencodeTag *ret);
/* Load a row from gencodeTag table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct gencodeTag *gencodeTagLoad(char **row);
/* Load a gencodeTag from row fetched with select * from gencodeTag
 * from database.  Dispose of this with gencodeTagFree(). */

struct gencodeTag *gencodeTagLoadAll(char *fileName);
/* Load all gencodeTag from whitespace-separated file.
 * Dispose of this with gencodeTagFreeList(). */

struct gencodeTag *gencodeTagLoadAllByChar(char *fileName, char chopper);
/* Load all gencodeTag from chopper separated file.
 * Dispose of this with gencodeTagFreeList(). */

#define gencodeTagLoadAllByTab(a) gencodeTagLoadAllByChar(a, '\t');
/* Load all gencodeTag from tab separated file.
 * Dispose of this with gencodeTagFreeList(). */

struct gencodeTag *gencodeTagCommaIn(char **pS, struct gencodeTag *ret);
/* Create a gencodeTag out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new gencodeTag */

void gencodeTagFree(struct gencodeTag **pEl);
/* Free a single dynamically allocated gencodeTag such as created
 * with gencodeTagLoad(). */

void gencodeTagFreeList(struct gencodeTag **pList);
/* Free a list of dynamically allocated gencodeTag's */

void gencodeTagOutput(struct gencodeTag *el, FILE *f, char sep, char lastSep);
/* Print out gencodeTag.  Separate fields with sep. Follow last field with lastSep. */

#define gencodeTagTabOut(el,f) gencodeTagOutput(el,f,'\t','\n');
/* Print out gencodeTag as a line in a tab-separated file. */

#define gencodeTagCommaOut(el,f) gencodeTagOutput(el,f,',',',');
/* Print out gencodeTag as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* GENCODETAG_H */

