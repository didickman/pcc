/*      $Id$    */
 /*
 * Copyright (c) 2020 Puresoftware Ltd.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * We define location operations which operate on the expression tree
 * during the first pass (before sending to the backend for code generation.)
 */

#include <assert.h>

#include "pass1.h"

#undef NIL
#define NIL NULL

#ifdef LANG_CXX
#define p1listf listf
#define p1tfree tfree
#else
#define NODE P1ND
#define talloc p1alloc
#define tcopy p1tcopy
#define nfree p1nfree
#undef n_type
#define n_type ptype
#undef n_qual
#define n_qual pqual
#undef n_df
#define n_df pdf
#define	n_ap pss
#endif

extern void defalign(int);

static char *
getsoname(struct symtab *sp)
{
	struct attr *ap;
	return (ap = attr_find(sp->sap, ATTR_SONAME)) ?
	    ap->sarg(0) : sp->sname;
}

#ifndef LANG_CXX
#define	sap sss
#endif

/*
 * clocal() is called to do local transformations on
 * an expression tree before being sent to the backend.
 */
NODE *
clocal(NODE *p)
{
	struct symtab *q;
	NODE *l, *r, *t;
	int o;
	int ty;
	int tmpnr, isptrvoid = 0;
	char *n;

	o = p->n_op;
	switch (o) {
		case STASG:

			l = p->n_left;
			r = p->n_right;
			if (r->n_op != STCALL && r->n_op != USTCALL)
				return p;

			/* assign left node as first argument to function */
			nfree(p);
			t = block(REG, NIL, NIL, r->n_type, r->n_df, r->n_ap);
			l->n_rval = R0;
			l = buildtree(ADDROF, l, NIL);
			l = buildtree(ASSIGN, t, l);

			if (r->n_right->n_op != CM) {
				r->n_right = block(CM, l, r->n_right, INT, 0, 0);
			} else {
				for (t = r->n_right; t->n_left->n_op == CM;
				    t = t->n_left);
				t->n_left = block(CM, l, t->n_left, INT, 0, 0);
			}
			return r;

		case CALL:
		case STCALL:
		case USTCALL:
			if (p->n_type == VOID)
				break;
			/*
			 * if the function returns void*, ecode() invokes
			 * delvoid() to convert it to uchar*.
			 * We just let this happen on the ASSIGN to the temp,
			 * and cast the pointer back to void* on access
			 * from the temp.
			 */
			if (p->n_type == PTR+VOID)
				isptrvoid = 1;
			r = tempnode(0, p->n_type, p->n_df, p->n_ap);
			tmpnr = regno(r);
			r = block(ASSIGN, r, p, p->n_type, p->n_df, p->n_ap);

			p = tempnode(tmpnr, r->n_type, r->n_df, r->n_ap);
			if (isptrvoid) {
				p = block(PCONV, p, NIL, PTR+VOID, p->n_df, 0);
			}
			p = buildtree(COMOP, r, p);
			break;

		case NAME:
			if ((q = p->n_sp) == NULL)
				return p;
			if (blevel == 0)
				return p;

			switch (q->sclass) {
				case PARAM:
				case AUTO:
					/* fake up a structure reference */
					r = block(REG, NIL, NIL, PTR+STRTY, 0, 0);
					slval(r, 0);
					r->n_rval = SPREG;
					p = stref(block(STREF, r, p, 0, 0, 0));
					break;
				case REGISTER:
					p->n_op = REG;
					slval(p, 0);
					p->n_rval = q->soffset;
					break;
				case STATIC:
					if (q->slevel > 0) {
						slval(p, 0);
						p->n_sp = q;
					}
					/* FALL-THROUGH */
				default:
					ty = p->n_type;
					n = getsoname(p->n_sp);
					if (strncmp(n, "__builtin", 9) == 0)
						break;
					p = block(ADDROF, p, NIL, INCREF(ty), p->n_df, p->n_ap);
					p = block(UMUL, p, NIL, ty, p->n_df, p->n_ap);
					break;
			}
			break;

		case STNAME:
			if ((q = p->n_sp) == NULL)
				return p;
			if (q->sclass != STNAME)
				return p;
			ty = p->n_type;
			p = block(ADDROF, p, NIL, INCREF(ty),p->n_df, p->n_ap);
			p = block(UMUL, p, NIL, ty, p->n_df, p->n_ap);
			break;

		case FORCE:
			/* put return value in return reg */
			p->n_op = ASSIGN;
			p->n_right = p->n_left;
			p->n_left = block(REG, NIL, NIL, p->n_type, 0, 0);
			p->n_left->n_rval = p->n_left->n_type == BOOL ? 
			    RETREG(BOOL_TYPE) : RETREG(p->n_type);
			break;

		case SCONV:
			l = p->n_left;
			if (p->n_type == l->n_type) {
				nfree(p);
				return l;
			}
			if ((p->n_type & TMASK) == 0 && (l->n_type & TMASK) == 0 &&
			    tsize(p->n_type, p->n_df, p->n_ap) == tsize(l->n_type, l->n_df, l->n_ap)) {
				if (p->n_type != FLOAT && p->n_type != DOUBLE &&
				    l->n_type != FLOAT && l->n_type != DOUBLE &&
				    l->n_type != LDOUBLE && p->n_type != LDOUBLE) {
					if (l->n_op == NAME || l->n_op == UMUL ||
					    l->n_op == TEMP) {
						l->n_type = p->n_type;
						nfree(p);
						return l;
					}
				}
			}

			if (l->n_op == ICON) {
				CONSZ val = glval(l);
				CONSZ lval;

				if (!ISPTR(p->n_type)) /* Pointers don't need to be conv'd */
				switch (p->n_type) {
					case BOOL:
						slval(l, glval(l) != 0);
						break;
					case CHAR:
						lval = (char)val;
						break;
					case UCHAR:
						lval = val & 0377;
					break;
				case SHORT:
					lval = (short)val;
					break;
				case USHORT:
					lval = val & 0177777;
					break;
				case ULONG:
				case UNSIGNED:
					lval = val & 0xffffffff;
					break;
				case LONG:
				case INT:
					lval = (int)val;
					break;
				case LONGLONG:
					lval = (long long)val;
					break;
				case ULONGLONG:
					lval = val;
					break;
				case VOID:
					break;
				case LDOUBLE:
				case DOUBLE:
				case FLOAT:
					break;
				default:
					cerror("unknown type %d", l->n_type);
				}
				if (p->n_type < FLOAT)
					slval(l, lval);
				l->n_type = p->n_type;
				l->n_ap = 0;
				nfree(p);
				return l;
			} else if (p->n_op == FCON) {
			cerror("SCONV FCON");
			}
			if ((DEUNSIGN(p->n_type) == CHAR ||
			    DEUNSIGN(p->n_type) == SHORT) &&
			    (l->n_type == FLOAT || l->n_type == DOUBLE ||
			    l->n_type == LDOUBLE)) {
				p = block(SCONV, p, NIL, p->n_type, p->n_df, p->n_ap);
				p->n_left->n_type = INT;
				return p;
			}
			break;

		case PCONV:
			l = p->n_left;
			if (l->n_op == ICON) {
				slval(l, (unsigned)glval(l));
				goto delp;
			}
			if (l->n_type < INT || DEUNSIGN(l->n_type) == LONGLONG) {
				p->n_left = block(SCONV, l, NIL, UNSIGNED, 0, 0);
				break;
			}
			if (l->n_op == SCONV)
				break;
			if (l->n_op == ADDROF && l->n_left->n_op == TEMP)
				goto delp;
			if (p->n_type > BTMASK && l->n_type > BTMASK)
				goto delp;
			break;

		delp:
			l->n_type = p->n_type;
			l->n_qual = p->n_qual;
			l->n_df = p->n_df;
			l->n_ap = p->n_ap;
			nfree(p);
			p = l;
			break;
	}

	return p;
}

/*
 * Called before sending the tree to the backend.
 */
void
myp2tree(NODE *p)
{
	struct symtab *sp;

	if (p->n_op != FCON)
		return;

#define IALLOC(sz)	(isinlining ? permalloc(sz) : tmpalloc(sz))

	sp = IALLOC(sizeof(struct symtab));
	sp->sclass = STATIC;
	sp->sap = 0;
	sp->slevel = 1; /* fake numeric label */
	sp->soffset = getlab();
	sp->sflags = 0;
	sp->stype = p->n_type;
	sp->squal = (CON >> TSHIFT);

	defloc(sp);
	inval(0, tsize(sp->stype, sp->sdf, sp->sap), p);

	p->n_op = NAME;
	slval(p, 0);
	p->n_sp = sp;
}

/*
 * Called during the first pass to determine if a NAME can be addressed.
 *
 * Return nonzero if supported, otherwise return 0.
 */
int
andable(NODE *p)
{
	if (blevel == 0)
		return 1;
	if (ISFTN(p->n_type))
		return 1;
	return 0;
}

/*
 * Return 1 if a variable of type 't' is OK to put in register.
 */
int
cisreg(TWORD t)
{
	if (t == FLOAT || t == DOUBLE || t == LDOUBLE)
		return 0; /* not yet */
	return 1;
}

/*
 * Allocate bits from the stack for dynamic-sized arrays.
 *
 * 'p' is the tree which represents the type being allocated.
 * 'off' is the number of 'p's to be allocated.
 * 't' is the storeable node where the address is written.
 */
void
spalloc(NODE *t, NODE *p, OFFSZ off)
{
	NODE *sp;

	p = buildtree(MUL, p, bcon(off/SZCHAR)); /* XXX word alignment? */

	/* sub the size from sp */
	sp = block(REG, NIL, NIL, p->n_type, 0, 0);
	slval(sp, 0);
	sp->n_rval = SP;
	ecomp(buildtree(MINUSEQ, sp, p));

	/* save the address of sp */
	sp = block(REG, NIL, NIL, PTR+INT, t->n_df, t->n_ap);
	slval(sp, 0);
	sp->n_rval = SP;
	t->n_type = sp->n_type;
	ecomp(buildtree(ASSIGN, t, sp));
}

/*
 * Print an integer constant node, may be associated with a label.
 * Do not free the node after use.
 * 'off' is bit offset from the beginning of the aggregate
 * 'fsz' is the number of bits this is referring to
 */
int
ninval(CONSZ off, int fsz, NODE *p)
{
	union { float f; double d; int i[2]; } u;
	struct symtab *q;
	TWORD t;
	int i, j;

	t = p->n_type;
	if (t > BTMASK)
		t = p->n_type = INT; /* pointer */

	if (p->n_op == ICON && p->n_sp != NULL && DEUNSIGN(t) != INT)
		uerror("element not constant");

	switch (t) {
		case LONGLONG:
		case ULONGLONG:
			i = (glval(p) >> 32);
			j = (glval(p) & 0xffffffff);
			p->n_type = INT;
			if (features(FEATURE_BIGENDIAN)) {
				slval(p, i);
				ninval(off+32, 32, p);
				slval(p, j);
				ninval(off, 32, p);
			} else {
				slval(p, j);
				ninval(off, 32, p);
				slval(p, i);
				ninval(off+32, 32, p);
			}
			break;
		case INT:
		case UNSIGNED:
			printf("\t.word 0x%x", (int)glval(p));
			if ((q = p->n_sp) != NULL) {
				if ((q->sclass == STATIC && q->slevel > 0)) {
					printf("+" LABFMT, q->soffset);
				} else
					printf("+%s", getexname(q));
			}
			printf("\n");
			break;
		default:
			return 0;
	}
	return 1;
}

/*
 * Prefix a leading underscore to a global variable (if necessary).
 */
char *
exname(char *p)
{
	return (p == NULL ? "" : p);
}

/*
 * Map types which are not defined on the local machine.
 */
TWORD
ctype(TWORD type)
{
	switch (BTYPE(type)) {
		case LONG:
			MODTYPE(type,LONGLONG);
			break;
		case ULONG:
			MODTYPE(type,ULONGLONG);
			break;
	}
	return (type);
}

/*
 * Before calling a function do any tree re-writing for the local machine.
 *
 * 'p' is the function tree (NAME)
 * 'q' is the CM-separated list of arguments.
 */
void
calldec(NODE *p, NODE *q) 
{
}

/*
 * While handling uninitialised variables, handle variables marked extern.
 */
void
extdec(struct symtab *q)
{
}

/* make a common declaration for id, if reasonable */
void
defzero(struct symtab *sp)
{
	int off;

	off = tsize(sp->stype, sp->sdf, sp->sap);
	off = (off+(SZCHAR-1))/SZCHAR;
	printf("	.%scomm ", sp->sclass == STATIC ? "l" : "");
	if (sp->slevel == 0)
		printf("%s,0%o\n", getexname(sp), off);
	else
		printf(LABFMT ",0%o\n", sp->soffset, off);
}

/*
 * va_start(ap, last) implementation.
 *
 * f is the NAME node for this builtin function.
 * a is the argument list containing:
 *	   CM
 *	ap   last
 */
NODE *
arm_builtin_stdarg_start(const struct bitable *bt, NODE *a)
{
	NODE *p, *q;
	int sz = 1;

	/* check num args and type */
	if (a == NULL || a->n_op != CM || a->n_left->n_op == CM ||
	    !ISPTR(a->n_left->n_type))
		goto bad;

	/* must first deal with argument size; use int size */
	p = a->n_right;
	if (p->n_type < INT) {
		/* round up to word */
		sz = SZINT / tsize(p->n_type, p->n_df, p->n_ap);
	}

	p = buildtree(ADDROF, p, NIL);  /* address of last arg */
	p = optim(buildtree(PLUS, p, bcon(sz)));
	q = block(NAME, NIL, NIL, PTR+VOID, 0, 0);
	q = buildtree(CAST, q, p);
	p = q->n_right;
	nfree(q->n_left);
	nfree(q);
	p = buildtree(ASSIGN, a->n_left, p);
	nfree(a);

	return p;

bad:
	uerror("bad argument to __builtin_stdarg_start");
	return bcon(0);
}

NODE *
arm_builtin_va_arg(const struct bitable *bt, NODE *a)
{
	NODE *p, *q, *r;
	int sz, tmpnr;

	/* check num args and type */
	if (a == NULL || a->n_op != CM || a->n_left->n_op == CM ||
	    !ISPTR(a->n_left->n_type) || a->n_right->n_op != TYPE)
		goto bad;

	r = a->n_right;

	/* get type size */
	sz = tsize(r->n_type, r->n_df, r->n_ap) / SZCHAR;
	if (sz < SZINT/SZCHAR) {
		werror("%s%s promoted to int when passed through ...",
			ISUNSIGNED(r->n_type) ? "unsigned " : "",
			DEUNSIGN(r->n_type) == SHORT ? "short" : "char");
		sz = SZINT/SZCHAR;
	}

	/* alignment */
	p = tcopy(a->n_left);
	if (sz > SZINT/SZCHAR && r->n_type != UNIONTY && r->n_type != STRTY) {
		p = buildtree(PLUS, p, bcon(ALSTACK/8 - 1));
		p = block(AND, p, bcon(-ALSTACK/8), p->n_type, p->n_df, p->n_ap);
	}

	/* create a copy to a temp node */
	q = tempnode(0, p->n_type, p->n_df, p->n_ap);
	tmpnr = regno(q);
	p = buildtree(ASSIGN, q, p);

	q = tempnode(tmpnr, p->n_type, p->n_df,p->n_ap);
	q = buildtree(PLUS, q, bcon(sz));
	q = buildtree(ASSIGN, a->n_left, q);

	q = buildtree(COMOP, p, q);

	nfree(a->n_right);
	nfree(a);

	p = tempnode(tmpnr, INCREF(r->n_type), r->n_df, r->n_ap);
	p = buildtree(UMUL, p, NIL);
	p = buildtree(COMOP, q, p);

	return p;

bad:
	uerror("bad argument to __builtin_va_arg");
	return bcon(0);
}

NODE *
arm_builtin_va_end(const struct bitable *bt, NODE *a)
{
	p1tfree(a);
 
	return bcon(0);
}

NODE *
arm_builtin_va_copy(const struct bitable *bt, NODE *a)
{
	NODE  *f;

	if (a == NULL || a->n_op != CM || a->n_left->n_op == CM)
		goto bad;
	f = buildtree(ASSIGN, a->n_left, a->n_right);
	nfree(a);
	return f;

bad:
	uerror("bad argument to __buildtin_va_copy");
	return bcon(0);
}

char *nextsect;
static int constructor;
static int destructor;

/*
 * Give target the opportunity of handling pragmas.
 */
int
mypragma(char *str)
{
	char *a2 = pragtok(NULL);

	if (strcmp(str, "tls") == 0) { 
		uerror("thread-local storage not supported for this target");
		return 1;
	} 
	if (strcmp(str, "constructor") == 0 || strcmp(str, "init") == 0) {
		constructor = 1;
		return 1;
	}
	if (strcmp(str, "destructor") == 0 || strcmp(str, "fini") == 0) {
		destructor = 1;
		return 1;
	}
	if (strcmp(str, "section") == 0 && a2 != NULL) {
		nextsect = newstring(a2, strlen(a2));
		return 1;
	}

	return 0;
}

/*
 * Called when a identifier has been declared, to give target last word.
 */
void
fixdef(struct symtab *sp)
{
	if ((constructor || destructor) && (sp->sclass != PARAM)) {
		printf("\t.section .%ctors,\"aw\",@progbits\n",
		    constructor ? 'c' : 'd');
		printf("\t.p2align 2\n");
		printf("\t.long %s\n", exname(sp->sname));
		printf("\t.previous\n");
		constructor = destructor = 0;
	}
}

void
pass1_lastchance(struct interpass *ip)
{
}
