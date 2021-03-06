/* altProbe.c was originally generated by the autoSql program, which also 
 * generated altProbe.h and altProbe.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "altProbe.h"


struct altProbe *altProbeLoad(char **row)
/* Load a altProbe from row fetched with select * from altProbe
 * from database.  Dispose of this with altProbeFree(). */
{
struct altProbe *ret;
int sizeOne;

AllocVar(ret);
ret->contProbeCount = sqlSigned(row[7]);
ret->alt1ProbeCount = sqlSigned(row[9]);
ret->alt2ProbeCount = sqlSigned(row[11]);
ret->transcriptCount = sqlSigned(row[13]);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlSigned(row[1]);
ret->chromEnd = sqlSigned(row[2]);
ret->type = sqlSigned(row[3]);
strcpy(ret->strand, row[4]);
ret->name = cloneString(row[5]);
ret->maxCounts = sqlSigned(row[6]);
sqlStringDynamicArray(row[8], &ret->contProbeSets, &sizeOne);
assert(sizeOne == ret->contProbeCount);
sqlStringDynamicArray(row[10], &ret->alt1ProbeSets, &sizeOne);
assert(sizeOne == ret->alt1ProbeCount);
sqlStringDynamicArray(row[12], &ret->alt2ProbeSets, &sizeOne);
assert(sizeOne == ret->alt2ProbeCount);
sqlStringDynamicArray(row[14], &ret->transcriptNames, &sizeOne);
assert(sizeOne == ret->transcriptCount);
return ret;
}

struct altProbe *altProbeLoadAll(char *fileName) 
/* Load all altProbe from a whitespace-separated file.
 * Dispose of this with altProbeFreeList(). */
{
struct altProbe *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[15];

while (lineFileNextCharRow(lf, '\t', row, 15))
    {
    el = altProbeLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct altProbe *altProbeLoadAllByChar(char *fileName, char chopper) 
/* Load all altProbe from a chopper separated file.
 * Dispose of this with altProbeFreeList(). */
{
struct altProbe *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[15];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = altProbeLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct altProbe *altProbeCommaIn(char **pS, struct altProbe *ret)
/* Create a altProbe out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new altProbe */
{
char *s = *pS;
int i;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlSignedComma(&s);
ret->chromEnd = sqlSignedComma(&s);
ret->type = sqlSignedComma(&s);
sqlFixedStringComma(&s, ret->strand, sizeof(ret->strand));
ret->name = sqlStringComma(&s);
ret->maxCounts = sqlSignedComma(&s);
ret->contProbeCount = sqlSignedComma(&s);
s = sqlEatChar(s, '{');
AllocArray(ret->contProbeSets, ret->contProbeCount);
for (i=0; i<ret->contProbeCount; ++i)
    {
    ret->contProbeSets[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
ret->alt1ProbeCount = sqlSignedComma(&s);
s = sqlEatChar(s, '{');
AllocArray(ret->alt1ProbeSets, ret->alt1ProbeCount);
for (i=0; i<ret->alt1ProbeCount; ++i)
    {
    ret->alt1ProbeSets[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
ret->alt2ProbeCount = sqlSignedComma(&s);
s = sqlEatChar(s, '{');
AllocArray(ret->alt2ProbeSets, ret->alt2ProbeCount);
for (i=0; i<ret->alt2ProbeCount; ++i)
    {
    ret->alt2ProbeSets[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
ret->transcriptCount = sqlSignedComma(&s);
s = sqlEatChar(s, '{');
AllocArray(ret->transcriptNames, ret->transcriptCount);
for (i=0; i<ret->transcriptCount; ++i)
    {
    ret->transcriptNames[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
*pS = s;
return ret;
}

void altProbeFree(struct altProbe **pEl)
/* Free a single dynamically allocated altProbe such as created
 * with altProbeLoad(). */
{
struct altProbe *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
/* All strings in contProbeSets are allocated at once, so only need to free first. */
if (el->contProbeSets != NULL)
    freeMem(el->contProbeSets[0]);
freeMem(el->contProbeSets);
/* All strings in alt1ProbeSets are allocated at once, so only need to free first. */
if (el->alt1ProbeSets != NULL)
    freeMem(el->alt1ProbeSets[0]);
freeMem(el->alt1ProbeSets);
/* All strings in alt2ProbeSets are allocated at once, so only need to free first. */
if (el->alt2ProbeSets != NULL)
    freeMem(el->alt2ProbeSets[0]);
freeMem(el->alt2ProbeSets);
/* All strings in transcriptNames are allocated at once, so only need to free first. */
if (el->transcriptNames != NULL)
    freeMem(el->transcriptNames[0]);
freeMem(el->transcriptNames);
freez(pEl);
}

void altProbeFreeList(struct altProbe **pList)
/* Free a list of dynamically allocated altProbe's */
{
struct altProbe *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    altProbeFree(&el);
    }
*pList = NULL;
}

void altProbeOutput(struct altProbe *el, FILE *f, char sep, char lastSep) 
/* Print out altProbe.  Separate fields with sep. Follow last field with lastSep. */
{
int i;
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%d", el->chromStart);
fputc(sep,f);
fprintf(f, "%d", el->chromEnd);
fputc(sep,f);
fprintf(f, "%d", el->type);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->strand);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%d", el->maxCounts);
fputc(sep,f);
fprintf(f, "%d", el->contProbeCount);
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->contProbeCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->contProbeSets[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
fputc(sep,f);
fprintf(f, "%d", el->alt1ProbeCount);
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->alt1ProbeCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->alt1ProbeSets[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
fputc(sep,f);
fprintf(f, "%d", el->alt2ProbeCount);
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->alt2ProbeCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->alt2ProbeSets[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
fputc(sep,f);
fprintf(f, "%d", el->transcriptCount);
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->transcriptCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->transcriptNames[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

