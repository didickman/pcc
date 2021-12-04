/*	$Id$	*/

/*
 * Copyright (c) 2021 Anders Magnusson. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <ctype.h>
#include <err.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>

#include "manifest.h"

/*
 * Peephole optimizer for Nova.
 * It is mandatory for this program to be run on the assembler file
 * after compilation, to cleanup short-distance jumps/loads/stores.
 * 
 */
#define BEGFTN	"#BEGFTN"
#define ENDFTN	"#ENDFTN"
#define BUFLEN	256

/*  constant lines */
struct symtab {
	DLIST_ENTRY(symtab) link;
	struct symtab *next;
	char *name;	/* e.g. "L345" */
	int off;	/* position in memory */
};

/* read instruction lines */
struct mblk {
	DLIST_ENTRY(mblk) link;
	char *name;	/* input string */
	int off;	/* insn offset in memory */
	int base;	/* long address position in memory (if needed) */
	struct symtab *sp;	/* symbol this jmp wants */
	int flags;
};

struct mblk mpole, mfirst, mlast;
struct symtab spole;

int lnum;
int curaddr;

static void readftn(void);
static void movecon(void);
static void printftn(void);

int
main(int argc, char *argv[])
{
	char buf[BUFLEN], *b;

	if (argc > 1)
		if (freopen(argv[1], "r", stdin) == NULL)
			err(1, "reopen stdin");
	if (argc > 2)
		if (freopen(argv[2], "w", stdout) == NULL)
			err(1, "reopen stdout");

	for (;;) {
		while ((b = fgets(buf, BUFLEN, stdin)) != NULL) {
			if (strncmp(b, BEGFTN, sizeof(BEGFTN) - 1) == 0)
				break;
			printf("%s", buf);
		}
		if (b == NULL)
			break;

		DLIST_INIT(&mpole,link);
		DLIST_INIT(&spole,link);
		DLIST_INIT(&mfirst,link);
		DLIST_INIT(&mlast,link);
		curaddr = 0;

		readftn();	/* fetch all function lines until end */
		movecon();	/* move constants to before/after functions */

		printftn();
	}
	return 0;
}

#define	IMCODE(a,b)	((a << 8) | b)
#define	IJMP	1
#define	IJSR	2
#define	IMREF	3
#define	IALC	4
static struct im {
	int code;
	int type;
} imatch[] = {
	{ IMCODE('j','m'), IJMP }, { IMCODE('j','s'), IJSR },
	{ IMCODE('l','d'), IMREF }, { IMCODE('s','t'), IMREF },
	{ IMCODE('i','s'), IMREF }, { IMCODE('d','s'), IMREF },
	{ IMCODE('m','o'), IALC },
	{ IMCODE('a','d'), IALC },
	{ IMCODE('s','u'), IALC },
	{ IMCODE('c','o'), IALC },
	{ IMCODE('n','e'), IALC },
	{ IMCODE('i','n'), IALC },
	{ IMCODE('a','n'), IALC },
	{ 0, 0 }
};

static int
itype(char *s)
{
	struct im *im = imatch;
	int c = (s[0] << 8) | s[1];

	for (; im->code; im++)
		if (im->code == c)
			return im->type;
	return 0;
}

static struct mblk *
mkmblk(char *s, int o)
{
	struct mblk *mb = malloc(sizeof(struct mblk));
	mb->name = s;
	mb->off = o;
	mb->base = mb->flags = 0;
	mb->sp = NULL;
	return mb;
}

static struct symtab *
mksp(char *s)
{
	struct symtab *sp = malloc(sizeof(struct symtab));
	sp->name = s;
	sp->off = 0;
	sp->next = 0;
	return sp;
}

/*
 * Readin function and store in linked list.
 * Check directives and do not save .globl.
 * This function expects the output syntax that pcc has, do not
 * try to parse anything else.
 */
static void
readftn(void)
{
	static struct symtab *savedsp;
	struct symtab *sp;
	struct mblk *mb;
	char buf[BUFLEN];
	char *b, *e;

	for (;;) {
		if ((b = fgets(buf, BUFLEN, stdin)) == NULL)
			err(1, "fgets");
		if (strncmp(b, ENDFTN, sizeof(ENDFTN)-1) == 0)
			return;
		while (*b == ' ' || *b == '\t')
			b++;
		for (e = b; *e; e++)
			;
		if (*--e != '\n')
			errx(1, "no lf line '%s'", buf);
		*e = 0; /* remove trailing lf */
		if (*b == '.') {
			if (strncmp(b, ".globl", 6) == 0) {
				printf("%s\n", buf);
				continue;
			} else if (strncmp(b, ".word", 5) != 0)
				errx(1, "bad dot dir %s", buf);
		}
		if (e[-1] == ':') {
			*--e = 0; /* remove : */
			sp = mksp(strdup(b));
			sp->next = savedsp;
			savedsp = sp;
			DLIST_INSERT_BEFORE(&spole, sp, link);
		} else {
			mb = mkmblk(strdup(b), curaddr++);
			mb->sp = savedsp;
			savedsp = NULL;
			DLIST_INSERT_BEFORE(&mpole, mb, link);
		}
	}
}

static struct mblk *
findw(char *s, struct mblk *pole)
{
	struct mblk *mb;

	DLIST_FOREACH(mb, pole, link) {
		if (strcmp(s, mb->name) == 0)
			return mb;
	}
	return NULL;
}

static int
isskip(char *s)
{
	while (*s && *s != ' ' && *s != '\t')
		s++;
	while (*s == ' ' || *s == '\t')
		s++;
	if (s[3] == ',')
		return 1;
	return 0;
}

/*
 * Move constants to a reasonable place.
 * The first 128 slots are put on mfirst list, the last 127 on mlast.
 * If any place inbetween is needed it is inserted on the correct place
 * in the mpole list.
 */
static void
movecon(void)
{
	struct mblk *mb, *mb2, *mblow, *mbhigh, mtemp, *mbpout;
	struct symtab *sp;
	char buf[BUFLEN];
	char *b, *s;
	int curcnt;

	if (curaddr > 128) {
		DLIST_FOREACH(mb, &mpole, link)
			if (mb->off == 128)
				break;
	} else
		mb = DLIST_PREV(&mpole, link);
	mblow = mb;
//printf("movecon1: mblow %p curaddr %d\n", mblow, curaddr);
	/*
	 * First handle constants at the beginning of the function.
	 * iterate backwards to get constants in the correct order.
	 */
	for (; mb != &mpole; mb = DLIST_PREV(mb, link)) {
		if ((b = strchr(mb->name, '['))) {
			s = strdup(b+1);
			s[strlen(s)-1] = 0; /* remove trailing ] */
			if ((mb2 = findw(s, &mfirst)) == 0) {
				mb2 = mkmblk(s, 0);
				DLIST_INSERT_AFTER(&mfirst, mb2, link);
			}
			if (mb2->sp == 0) {
				sprintf(buf, "LP%d", lnum++);
				sp = mksp(strdup(buf));
				mb2->sp = sp;
				DLIST_INSERT_BEFORE(&spole, sp, link);
			}
			strcpy(b, mb2->sp->name);
		}
	}
	if (curaddr <= 128)
		return; /* <= 128 entries */
//printf("movecon2: mblow %p curaddr %d\n", mblow, curaddr);
	/*
	 * Second; move constants at the end of the function.
	 */
	if (curaddr <= 255) {
		mb = DLIST_NEXT(mblow, link);
	} else {
		DLIST_FOREACH_REVERSE(mb, &mpole, link)
			if (mb->off == curaddr-127)
				break;
	}
	mbhigh = mb;
//printf("movecon3: mbhigh %p curaddr %d\n", mbhigh, curaddr);
	for (; mb != &mpole; mb = DLIST_NEXT(mb, link)) {
		if ((b = strchr(mb->name, '['))) {
			s = strdup(b+1);
			s[strlen(s)-1] = 0; /* remove trailing ] */
			if ((mb2 = findw(s, &mlast)) == 0) {
				mb2 = mkmblk(s, 0);
				DLIST_INSERT_BEFORE(&mlast, mb2, link);
			}
			if (mb2->sp == 0) {
				sprintf(buf, "LP%d", lnum++);
				sp = mksp(strdup(buf));
				mb2->sp = sp;
				DLIST_INSERT_BEFORE(&spole, sp, link);
			}
			strcpy(b, mb2->sp->name);
		}
	}
//printf("movecon5: mb %p curaddr %d\n", mb, curaddr);
	if (curaddr <= 255)
		return;
//printf("movecon6: mb %p curaddr %d\n", mb, curaddr);
	/*
	 * Here begins the tricky part.
	 * - Start at mblow. At the first memref insn start counting.
	 * - If a valid jmp is found, write out the collected data.
	 * - if passed 127 insns (or reached mbhigh) and no valid jmp,
	 *   go back (at least) one instruction, insert a jmp,
	 *   and write out the collected data. if passed mbhigh return.
	 * - go back and try again.
	 */
	DLIST_INIT(&mtemp, link);
	curcnt = 0;
	mbpout = NULL;

//printf("movecon7: mb %p curaddr %d\n", mb, curaddr);
	for (mb = DLIST_NEXT(mblow, link); mb != mbhigh;
	    mb = DLIST_NEXT(mb, link)) {
		int ityp = itype(mb->name);

//printf("movecon4: mb %p mbnum %d lnum %d line %s\n", mb, mb->off, curcnt, mb->name);
		if (curcnt)
			curcnt++;
		if (ityp == IMREF && (b = strchr(mb->name, '['))) {
			s = strdup(b+1);
			s[strlen(s)-1] = 0;	/* remove trailing ] */
			if ((mb2 = findw(s, &mtemp)) == 0) {
				mb2 = mkmblk(s, 0);
				DLIST_INSERT_AFTER(&mtemp, mb2, link);
			}
			if (mb2->sp == 0) {
				sprintf(buf, "LP%d", lnum++);
				sp = mksp(strdup(buf));
				mb2->sp = sp;
				DLIST_INSERT_BEFORE(&spole, sp, link);
			}
			strcpy(b, mb2->sp->name);
			if (curcnt == 0)
				curcnt++;
		}
		if (curcnt && ityp == IJMP) {
			mb2 = DLIST_PREV(mb, link);
			if (itype(mb2->name) != IALC || !isskip(mb2->name)) {
				mbpout = DLIST_NEXT(mb, link);
				/* safe to print out words after jmp */
				while (!DLIST_ISEMPTY(&mtemp, link)) {
					mb2 = DLIST_NEXT(&mtemp, link);
					DLIST_REMOVE(mb2, link);
					DLIST_INSERT_AFTER(mb, mb2, link);
				}
				curcnt = 0;
			}
		}
		if (curcnt == 127)
			errx(1, "curcnt == 127");
	}
	if (curcnt) {
		if (mb->off - mbpout->off < 128) {
			/* can put constants at last position */
			while (!DLIST_ISEMPTY(&mtemp, link)) {
				mb2 = DLIST_PREV(&mtemp, link);
				DLIST_REMOVE(mb2, link);
				DLIST_INSERT_BEFORE(mbpout, mb2, link);
			}
		} else
			errx(1, "FIXME long distance");
	}
//printf("movecon8: mb %p mbnum %d\n", mb, mb->off);
}

/*
 * When everything is modified; print out the new function.
 */
static void
printftn(void)
{
	struct symtab *sp;
	struct mblk *mb;

	DLIST_FOREACH(mb, &mfirst, link)
		printf("%s:\t%s\n", mb->sp->name, mb->name);

	DLIST_FOREACH(mb, &mpole, link) {
		while ((sp = mb->sp) != NULL) {
			printf("%s:\n", sp->name);
			mb->sp = sp->next;
		}
		printf("	%s\n", mb->name);
	}

	DLIST_FOREACH(mb, &mlast, link)
		printf("%s:\t%s\n", mb->sp->name, mb->name);

}
