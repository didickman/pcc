/*	$Id$	*/
/*
 * Copyright (c) 2003 Anders Magnusson (ragge@ludd.luth.se).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission
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


# include "pass2.h"

# define TLL TLONGLONG|TULONGLONG
# define ANYSIGNED TINT|TLONG|TSHORT|TCHAR
# define ANYUSIGNED TUNSIGNED|TULONG|TUSHORT|TUCHAR
# define ANYFIXED ANYSIGNED|ANYUSIGNED
# define TUWORD TUNSIGNED
# define TSWORD TINT
# define TWORD TUWORD|TSWORD
# define ANYSH	SCON|SAREG|SOREG|SNAME
# define ARONS	SAREG|SOREG|SNAME|STARNM

#define XSL(c)  NEEDS(NREG(c, 1), NSL(c))
#define NAREG   NEEDS(NREG(A, 1))
#define NBREG   NEEDS(NREG(B, 1))
#define NCREG   NEEDS(NREG(C, 1))   

struct optab table[] = {
/* First entry must be an empty entry */
{ -1, FOREFF, SANY, TANY, SANY, TANY, 0, 0, "", },

/* PCONVs are usually not necessary */
{ PCONV,	INAREG,
	SAREG,	TWORD|TPOINT,
	SAREG,	TWORD|TPOINT,
		0,	RLEFT,
		"", },

/* long -> ptr */
{ PCONV,	INAREG,
	SBREG|SOREG|SNAME,	TLONG|TULONG,
	SAREG,			TPOINT,
		XSL(A),	RESC1,
		"mov	UL,A1\n", },

/* convert uchar to char; sign-extend byte */
{ SCONV,	INAREG,
	SAREG,	TUCHAR,
	SAREG,	TCHAR,
		0,	RLEFT,
		"movb	AL,AL\n", },

/* convert char to uchar; zero-extend byte */
{ SCONV,	INAREG,
	SAREG,	TCHAR,
	SAREG,	TUCHAR,
		0,	RLEFT,
		"bic	$177400,AL\n", },

/* convert char to int or unsigned.  Already so in regs */
{ SCONV,	INAREG,
	SAREG,	TCHAR,
	SAREG,	TINT|TUNSIGNED,
		0,	RLEFT,
		"", },

/* char (in mem) to (u)int */
{ SCONV,	INAREG,
	SOREG|SCON|SNAME,	TCHAR,
	SAREG,	TINT,
		XSL(A),	RESC1,
		"movb	AL,A1\n", },

/* convert uchar to (u)int in reg */
/* Nothing to do since reg is already large enough */
{ SCONV,	INAREG,
	SAREG,	TUCHAR,
	SAREG,	TINT|TUNSIGNED,
		0,	RLEFT,
		"", },

/* convert uchar to (u)int from mem */
{ SCONV,	INAREG,
	SOREG|SNAME,	TUCHAR,
	SAREG,	TINT|TUNSIGNED,
		NAREG,	RESC1,
		"clr	A1\nbisb	AL,A1\n", },

/* convert char to (u)long */
{ SCONV,	INBREG,
	SAREG|SOREG|SCON|SNAME,	TCHAR,
	SANY,	TULONG|TLONG,
		NBREG,	RESC1,
		"movb	AL,U1\nsxt	A1\n", },

/* convert uchar to ulong */
{ SCONV,	INBREG,
	SAREG,	TUCHAR,
	SANY,	TLONG|TULONG,
		XSL(B),	RESC1,
		"mov	AL,U1\nclr	A1\n", },

/* (u)char -> float/double */
{ SCONV,	INCREG,
	SAREG,	TCHAR|TUCHAR,
	SANY,	TFLOAT|TDOUBLE,
		NCREG,	RESC1,
		"movif	AL,A1\n", },

/* convert (u)int to char */
{ SCONV,	INAREG,
	SAREG,	TWORD,
	SANY,	TCHAR,
		0,	RLEFT,
		"movb	AL,AL\n", },

/* convert (u)int to uchar  */
{ SCONV,	INAREG,
	SAREG,	TWORD,
	SANY,	TUCHAR,
		NAREG,	RESC1,
		"clr	A1\nbisb	AL,A1\n", },

/* convert (u)int to (u)int */
{ SCONV,	INAREG,
	SAREG,	TWORD,
	SANY,	TWORD,
		0,	RLEFT,
		"", },

/* convert pointer to char */
{ SCONV,	INAREG,
	SAREG,	TPOINT,
	SANY,	TCHAR,
		0,	RLEFT,
		"movb	AL,AL\n", },

/* convert pointer to (u)int */
{ SCONV,	INAREG,
	SAREG,	TPOINT,
	SANY,	TWORD,
		0,	RLEFT,
		"", },

/* convert int to long from memory */
{ SCONV,	INBREG,
	SNAME|SOREG,	TINT,
	SANY,	TLONG,
		NBREG,	RESC1,
		"mov	AL,U1\nsxt	A1\n", },

/* int -> (u)long. XXX - only in r0 and r1 */
{ SCONV,	INBREG,
	SAREG,	TINT|TPOINT,
	SANY,	TLONG|TULONG,
		NEEDS(NREG(B, 1), NSL(B), NLEFT(R1), NRES(R01)), RESC1,
		"tst	AL\nsxt	r0\n", },

/* int -> float/double */
{ SCONV,	INCREG,
	SAREG|SNAME|SOREG,	TINT,
	SANY,	TFLOAT|TDOUBLE,
		NCREG,	RESC1,
		"movif	AL,A1\n", },


/* unsigned -> (u)long. XXX - only in r0 and r1 */
{ SCONV,	INBREG,
	SAREG,	TUNSIGNED,
	SANY,	TLONG|TULONG,
		NEEDS(NREG(B, 1), NSL(B), NLEFT(R1), NRES(R01)), RESC1,
		"clr	r0\n", },

/* uint -> double */
{ SCONV,	INCREG,
	SAREG|SNAME|SOREG|SCON,	TUNSIGNED,
	SANY,			TFLOAT|TDOUBLE,
		XSL(C),	RESC1,
		"mov	AL,-(sp)\nclr	-(sp)\n"
		"setl\nmovif	(sp)+,A1\nseti\n", },

/* (u)long -> char */
{ SCONV,	INAREG,
	SBREG|SOREG|SNAME,	TLONG|TULONG,
	SAREG,			TCHAR,
		XSL(A),	RESC1,
		"movb	UL,A1\n", },

/* (u)long -> uchar */
{ SCONV,	INAREG,
	SBREG|SOREG|SNAME,	TLONG|TULONG,
	SAREG,			TUCHAR,
		XSL(A),	RESC1,
		"clr	A1\nbisb	UL,A1\n", },

/* long -> int */
{ SCONV,	INAREG,
	SBREG|SOREG|SNAME,	TLONG|TULONG,
	SAREG,			TWORD,
		XSL(A),	RESC1,
		"mov	UL,A1\n", },

/* (u)long -> (u)long, nothing */
{ SCONV,	INBREG,
	SBREG,	TLONG|TULONG,
	SANY,	TLONG|TULONG,
		XSL(B),	RESC1,
		"", },

/* long -> float/double */
{ SCONV,	INCREG,
	SBREG|SNAME|SOREG|SCON,	TLONG,
	SANY,		TFLOAT|TDOUBLE,
		XSL(C),	RESC1,
		"mov	UL,-(sp)\nmov	AL,-(sp)\n"
		"setl\nmovif	(sp)+,A1\nseti\n", },

/* ulong -> float/double */
{ SCONV,	INCREG,
	SBREG|SNAME|SOREG|SCON,	TULONG,
	SANY,		TFLOAT|TDOUBLE,
		XSL(C),	RESC1,
		"mov	UL,(sp)\nmov	AL,-(sp)\n"
		"jsr	pc,ultof\ntst	(sp)+\n", },

/* float/double -> (u)long XXX - correct? */
{ SCONV,	INBREG,
	SCREG,	TFLOAT|TDOUBLE,
	SANY,	TLONG|TULONG,
		NBREG,	RESC1,
		"setl\nmovfi	AL,-(sp)\nmov	(sp)+,A1\n"
		"mov	(sp)+,U1\nseti\n", },

/* double -> int*/
{ SCONV,	INAREG,
	SCREG,	TFLOAT|TDOUBLE,
	SANY,	TINT,
		NAREG,	RESC1,
		"movfi	AL,A1\n", },

/* float/double -> float/double */
/* In regs. floating point always stored as double in regs */
{ SCONV,	INCREG,
	SCREG,	TFLOAT|TDOUBLE,
	SANY,	TANY,
		0,	RLEFT,
		"", },

/*
 * Subroutine calls.
 */
{ CALL,		INAREG,
	SCON,	TANY,
	SAREG,	TWORD|TPOINT|TCHAR|TUCHAR,
		XSL(A),	RESC1,
		"jsr	pc,*AL\nZC", },

{ UCALL,	INAREG,
	SCON,	TANY,
	SAREG,	TWORD|TPOINT|TCHAR|TUCHAR,
		XSL(A),	RESC1,
		"jsr	pc,*AL\n", },

{ CALL,		INAREG,
	SAREG,	TANY,
	SANY,	TANY,
		XSL(A),	RESC1,	/* should be 0 */
		"jsr	pc,(AL)\nZC", },

{ UCALL,	INAREG,
	SAREG,	TANY,
	SANY,	TANY,
		XSL(A),	RESC1,	/* should be 0 */
		"jsr	pc,(AL)\n", },

{ CALL,		INBREG,
	SCON,	TANY,
	SBREG,	TLONG|TULONG,
		XSL(B),	RESC1,
		"jsr	pc,*AL\nZC", },

{ UCALL,	INBREG,
	SCON,	TANY,
	SBREG,	TLONG|TULONG,
		XSL(B),	RESC1,
		"jsr	pc,*AL\n", },

{ CALL,		INBREG,
	SAREG,	TANY,
	SBREG,	TLONG|TULONG,
		XSL(B),	RESC1,
		"jsr	pc,(AL)\nZC", },

{ UCALL,	INBREG,
	SAREG,	TANY,
	SBREG,	TLONG|TULONG,
		XSL(B),	RESC1,
		"jsr	pc,(AL)\n", },

{ CALL,		INCREG,
	SCON,	TANY,
	SCREG,	TFLOAT|TDOUBLE,
		NCREG,	RESC1,
		"jsr	pc,*AL\nZC", },

{ UCALL,	INCREG,
	SCON,	TANY,
	SCREG,	TFLOAT|TDOUBLE,
		NCREG,	RESC1,
		"jsr	pc,*AL\n", },

{ CALL,		INCREG,
	SAREG,	TANY,
	SCREG,	TFLOAT|TDOUBLE,
		XSL(C),	RESC1,
		"jsr	pc,(AL)\nZC", },

{ UCALL,	INCREG,
	SAREG,	TANY,
	SCREG,	TFLOAT|TDOUBLE,
		XSL(C),	RESC1,
		"jsr	pc,(AL)\n", },

{ CALL,		FOREFF,
	SCON,	TANY,
	SANY,	TANY,
		0,	0,
		"jsr	pc,*AL\nZC", },

{ UCALL,	FOREFF,
	SCON,	TANY,
	SANY,	TANY,
		0,	0,
		"jsr	pc,*AL\n", },

{ CALL,		FOREFF,
	SAREG,	TANY,
	SANY,	TANY,
		0,	0,
		"jsr	pc,(AL)\nZC", },

{ UCALL,	FOREFF,
	SAREG,	TANY,
	SANY,	TANY,
		0,	0,
		"jsr	pc,(AL)\n", },

{ STCALL,	INAREG,
	SCON|SOREG|SNAME,	TANY,
	SANY,	TANY,
		XSL(A),	RESC1,
		"jsr	pc,*AL\nZC", },

{ USTCALL,	INAREG,
	SCON|SOREG|SNAME,	TANY,
	SANY,	TANY,
		XSL(A),	RESC1,
		"jsr	pc,*AL\n", },

{ STCALL,	FOREFF,
	SCON|SOREG|SNAME,	TANY,
	SANY,	TANY,
		0,	0,
		"jsr	pc,*AL\nZC", },

{ USTCALL,	FOREFF,
	SCON|SOREG|SNAME,	TANY,
	SANY,	TANY,
		0,	0,
		"jsr	pc,*AL\n", },

{ STCALL,	INAREG,
	SAREG,	TANY,
	SANY,	TANY,
		XSL(A),	RESC1,	/* should be 0 */
		"jsr	pc,(AL)\nZC", },

{ USTCALL,	INAREG,
	SAREG,	TANY,
	SANY,	TANY,
		XSL(A),	RESC1,	/* should be 0 */
		"jsr	pc,(AL)\n", },

/*
 * The next rules handle all binop-style operators.
 */
/* Add one to anything left but use only for side effects */
{ PLUS,		FOREFF|INAREG|FORCC,
	SAREG|SNAME|SOREG,	TWORD|TPOINT,
	SONE,			TANY,
		0,	RLEFT|RESCC,
		"inc	AL\n", },

/* add one for char to reg, special handling */
{ PLUS,		FOREFF|INAREG|FORCC,
	SAREG,	TCHAR|TUCHAR,
	SONE,		TANY,
		0,	RLEFT|RESCC,
		"inc	AL\n", },

/* add one for char to memory */
{ PLUS,		FOREFF|FORCC,
	SNAME|SOREG|STARNM,	TCHAR|TUCHAR,
	SONE,			TANY,
		0,	RLEFT|RESCC,
		"incb	AL\n", },

{ PLUS,		INBREG|FOREFF,
	SBREG,			TLONG|TULONG,
	SBREG|SNAME|SOREG|SCON,	TLONG|TULONG,
		0,	RLEFT,
		"add	AR,AL\nadd	UR,UL\nadc	AL\n", },

/* Integer to pointer addition */
{ PLUS,		INAREG,
	SAREG,	TPOINT|TWORD,
	SAREG,	TINT|TUNSIGNED,
		0,	RLEFT,
		"add	AR,AL\n", },

/* Add to reg left and reclaim reg */
{ PLUS,		INAREG|FOREFF|FORCC,
	SAREG|SNAME|SOREG,	TWORD|TPOINT,
	SAREG|SNAME|SOREG|SCON,	TWORD|TPOINT,
		0,	RLEFT|RESCC,
		"add	AR,AL\n", },

/* Add to anything left but use only for side effects */
{ PLUS,		FOREFF|FORCC,
	SNAME|SOREG,	TWORD|TPOINT,
	SAREG|SNAME|SOREG|SCON,	TWORD|TPOINT,
		0,	RLEFT|RESCC,
		"add	AR,AL\n", },

{ PLUS,		INAREG|FOREFF|FORCC,
	SAREG,			TCHAR|TUCHAR,
	SAREG|SNAME|SOREG|SCON,	TCHAR|TUCHAR,
		0,	RLEFT|RESCC,
		"add	AR,AL\n", },

/* floating point */
{ PLUS,		INCREG|FOREFF|FORCC,
	SCREG,	TFLOAT,
	SCREG,	TFLOAT,
		0,	RLEFT|RESCC,
		"addf	AR,AL\n", },

{ PLUS,		INCREG|FOREFF|FORCC,
	SCREG,			TDOUBLE,
	SCREG|SNAME|SOREG,	TDOUBLE,
		0,	RLEFT|RESCC,
		"addf	AR,AL\n", },

/* Post-increment read, byte */
{ MINUS,	INAREG,
	SINCB,	TCHAR|TUCHAR,
	SONE,	TANY,
		NAREG,	RESC1,
		"movb	ZG,A1\nincb	ZG\n", },

/* Post-increment read, int */
{ MINUS,	INAREG,
	SINCB,	TWORD|TPOINT,
	SONE,	TANY,
		NAREG,	RESC1,
		"mov	ZG,A1\ninc	ZG\n", },

{ MINUS,		INBREG|FOREFF,
	SBREG,			TLONG|TULONG,
	SBREG|SNAME|SOREG|SCON,	TLONG|TULONG,
		0,	RLEFT,
		"sub	AR,AL\nsub	UR,UL\nsbc	AL\n", },

/* Sub one from anything left */
{ MINUS,	FOREFF|INAREG|FORCC,
	SAREG|SNAME|SOREG,	TWORD|TPOINT,
	SONE,			TANY,
		0,	RLEFT|RESCC,
		"dec	AL\n", },

{ MINUS,		INAREG|FOREFF,
	SAREG,			TWORD|TPOINT,
	SAREG|SNAME|SOREG|SCON,	TWORD|TPOINT,
		0,	RLEFT,
		"sub	AR,AL\n", },

/* Sub from anything left but use only for side effects */
{ MINUS,	FOREFF|INAREG|FORCC,
	SAREG|SNAME|SOREG,	TWORD|TPOINT,
	SAREG|SNAME|SOREG|SCON,	TWORD|TPOINT,
		0,	RLEFT|RESCC,
		"sub	AR,AL\n", },

/* Sub one left but use only for side effects */
{ MINUS,	FOREFF|FORCC,
	SAREG|SNAME|SOREG,	TCHAR|TUCHAR,
	SONE,			TANY,
		0,	RLEFT|RESCC,
		"decb	AL\n", },

/* Sub from anything left but use only for side effects */
{ MINUS,		FOREFF|FORCC,
	SAREG|SNAME|SOREG,	TCHAR|TUCHAR,
	SAREG|SNAME|SOREG|SCON,	TCHAR|TUCHAR|TWORD|TPOINT,
		0,	RLEFT|RESCC,
		"subb	AR,AL\n", },

/* floating point */
{ MINUS,	INCREG|FOREFF|FORCC,
	SCREG,	TFLOAT,
	SCREG,	TFLOAT,
		0,	RLEFT|RESCC,
		"subf	AR,AL\n", },

{ MINUS,	INCREG|FOREFF|FORCC,
	SCREG,			TDOUBLE,
	SCREG|SNAME|SOREG,	TDOUBLE,
		0,	RLEFT|RESCC,
		"subf	AR,AL\n", },

/*
 * The next rules handle all shift operators.
 */
{ LS,	INBREG|FOREFF,
	SBREG,	TLONG|TULONG,
	SAREG,	TINT|TUNSIGNED,
		0,	RLEFT,
		"ashc	AR,AL\n", },

{ LS,	INAREG|FOREFF,
	SAREG,	TWORD,
	SONE,	TANY,
		0,	RLEFT,
		"asl	AL\n", },

{ LS,	INAREG|FOREFF,
	SAREG,	TINT|TCHAR|TUNSIGNED|TUCHAR,
	ANYSH,	TWORD,
		0,	RLEFT,
		"ash	AR,AL\n", },

{ RS,	INAREG|FOREFF,
	SAREG,	TUNSIGNED|TUCHAR,
	SAREG|SCON,	TWORD,
		NEEDS(NEVER(R0), NEVER(R1), NLEFT(R1), NRES(R1), 
		    NORIGHT(R0), NORIGHT(R1)),	RLEFT,
		"clr	r0\nashc	AR,r0\n", },

{ RS,	INAREG|FOREFF,
	SAREG,	TINT|TCHAR,
	SAREG|SCON,	TWORD,
		0,	RLEFT,
		"ash	AR,AL\n", },

/*
 * The next rules takes care of assignments. "=".
 */

/* First optimizations, in lack of weight it uses first found */
/* Start with class A registers */

/* Clear word at address */
{ ASSIGN,	FOREFF|FORCC,
	ARONS,	TWORD|TPOINT,
	SZERO,		TANY,
		0,	RESCC,
		"clr	AL\n", },

/* Clear word at reg */
{ ASSIGN,	FOREFF|INAREG,
	SAREG,	TWORD|TPOINT,
	SZERO,		TANY,
		0,	RDEST,
		"clr	AL\n", },

/* Clear byte at address.  No reg here. */
{ ASSIGN,	FOREFF,
	SNAME|SOREG|STARNM,	TCHAR|TUCHAR,
	SZERO,		TANY,
		0,	RDEST,
		"clrb	AL\n", },

/* Clear byte in reg. must clear the whole register. */
{ ASSIGN,	FOREFF|INAREG,
	SAREG,	TCHAR|TUCHAR,
	SZERO,	TANY,
		0,	RDEST,
		"clr	AL\n", },

/* The next is class B regs */

/* Clear long at address or in reg */
{ ASSIGN,	FOREFF|INBREG,
	SNAME|SOREG|SBREG,	TLONG|TULONG,
	SZERO,			TANY,
		0,	RDEST,
		"clr	AL\nclr	UL\n", },

/* Save 2 bytes if high-order bits are zero */
{ ASSIGN,	FOREFF|INBREG,
	SBREG,	TLONG|TULONG,
	SSCON,	TLONG,
		0,	RDEST,
		"mov	UR,UL\nsxt	AL\n", },

/* Must have multiple rules for long otherwise regs may be trashed */
{ ASSIGN,	FOREFF|INBREG,
	SBREG,			TLONG|TULONG,
	SCON|SNAME|SOREG,	TLONG|TULONG,
		0,	RDEST,
		"mov	AR,AL\nmov	UR,UL\n", },

{ ASSIGN,	FOREFF|INBREG,
	SNAME|SOREG,	TLONG|TULONG,
	SBREG,			TLONG|TULONG,
		0,	RDEST,
		"mov	AR,AL\nmov	UR,UL\n", },

{ ASSIGN,	FOREFF,
	SNAME|SOREG,	TLONG|TULONG,
	SCON|SNAME|SOREG,	TLONG|TULONG,
		0,	0,
		"mov	AR,AL\nmov	UR,UL\n", },

{ ASSIGN,	INBREG|FOREFF,
	SBREG,	TLONG|TULONG,
	SBREG,	TLONG|TULONG,
		0,	RDEST,
		"ZE\n", },

{ ASSIGN,	FOREFF|INAREG|FORCC,
	SAREG,			TWORD|TPOINT,
	SAREG|SNAME|SOREG|SCON,	TWORD|TPOINT,
		0,	RDEST|RESCC,
		"mov	AR,AL\n", },

{ ASSIGN,	FOREFF|INAREG|FORCC,
	ARONS,	TWORD|TPOINT,
	SAREG,	TWORD|TPOINT,
		0,	RDEST|RESCC,
		"mov	AR,AL\n", },

{ ASSIGN,	FOREFF|FORCC,
	SNAME|SOREG,		TWORD|TPOINT,
	SNAME|SOREG|SCON,	TWORD|TPOINT,
		0,	RESCC,
		"mov	AR,AL\n", },

{ ASSIGN,	FOREFF|INAREG|FORCC,
	SAREG,		TCHAR|TUCHAR,
	ARONS|SCON,	TCHAR|TUCHAR,
		0,	RDEST|RESCC,
		"movb	AR,AL\n", },

{ ASSIGN,	FOREFF|INAREG|FORCC,
	ARONS,	TCHAR|TUCHAR,
	SAREG,	TCHAR|TUCHAR,
		0,	RDEST|RESCC,
		"movb	AR,AL\n", },

{ ASSIGN,	FOREFF|FORCC,
	SNAME|SOREG|STARNM,		TCHAR|TUCHAR,
	SNAME|SOREG|SCON|STARNM,	TCHAR|TUCHAR,
		0,	RDEST|RESCC,
		"movb	AR,AL\n", },

/* Floating point */
{ ASSIGN,	FOREFF|INCREG,
	SCREG,		TDOUBLE,
	SNAME|SOREG|SCON,	TDOUBLE,
		0,	RDEST,
		"movf	AR,AL\n", },

{ ASSIGN,	FOREFF|INCREG,
	SCREG,		TFLOAT,
	SNAME|SOREG|SCON,	TFLOAT,
		0,	RDEST,
		"movof	AR,AL\n", },

{ ASSIGN,	FOREFF|INCREG,
	SNAME|SOREG|SCREG,	TDOUBLE,
	SCREG,			TDOUBLE,
		0,	RDEST,
		"movf	AR,AL\n", },

{ ASSIGN,	FOREFF|INCREG,
	SNAME|SOREG|SCREG,	TFLOAT,
	SCREG,			TFLOAT,
		0,	RDEST,
		"movfo	AR,AL\n", },

/* Struct assigns */
{ STASG,	FOREFF|INAREG,
	SOREG|SNAME,	TANY,
	SAREG,		TPTRTO|TANY,
		NEEDS(NEVER(R0), NEVER(R1), NRIGHT(R0), NOLEFT(R0),
		    NOLEFT(R2), NORIGHT(R2), NEVER(R2)),	RDEST,
		"Fpush	r0\nZIFpop	r0\n", },

/*
 * DIV/MOD/MUL 
 */
/* XXX - mul may use any odd register, only r1 for now */
{ MUL,	INAREG,
	SAREG,				TWORD|TPOINT,
	SAREG|SNAME|SOREG|SCON,		TWORD|TPOINT,
		NEEDS(NLEFT(R1), NRES(R1)),	RLEFT,
		"mul	AR,AL\n", },

{ MUL,	INBREG,
	SBREG|SNAME|SCON|SOREG,		TLONG|TULONG,
	SBREG|SNAME|SCON|SOREG,		TLONG|TULONG,
		NEEDS(NREG(B, 1), NSL(B), NSR(B), NEVER(R0),
		    NEVER(R1), NRES(R01)), 	RESC1,
#if 0
		"mov	UR,-(sp)\nmov	AR,-(sp)\n"
		"mov	UL,-(sp)\nmov	AL,-(sp)\n"
		"jsr	pc,lmul\nadd	$10,sp\n",
#endif
		"ZK", },

{ MUL,	INCREG,
	SCREG,		TFLOAT,
	SCREG,		TFLOAT,
		0,	RLEFT,
		"mulf	AR,AL\n", },

{ MUL,	INCREG,
	SCREG,			TDOUBLE,
	SCREG|SNAME|SOREG,	TDOUBLE,
		0,	RLEFT,
		"mulf	AR,AL\n", },

/* need extra move to be sure N flag is correct for sxt */
{ DIV,	INAREG,
	SAREG,				TINT|TPOINT,
	SAREG|SNAME|SCON|SOREG,		TINT|TPOINT,
		NEEDS(NEVER(R0), NEVER(R1), NLEFT(R1),
		    NRES(R0), NORIGHT(R0), NORIGHT(R1)), 	RDEST,
		"tst	r1\nsxt	r0\ndiv	AR,r0\n", },

/* udiv uses args in registers */
{ DIV,	INAREG,
	SAREG,		TUNSIGNED,
	SAREG,		TUNSIGNED,
		NEEDS(NREG(A, 1), NSL(A), NSR(A), NLEFT(R0), 
		    NRIGHT(R1), NRES(R0)),	RESC1,
		"jsr	pc,udiv\n", },

{ DIV,	INBREG,
	SBREG|SNAME|SCON|SOREG,		TLONG|TULONG,
	SBREG|SNAME|SCON|SOREG,		TLONG|TULONG,
		NEEDS(NREG(B, 1), NSL(B), NSR(B), NEVER(R0), 
		    NEVER(R1), NRES(R01)),	RESC1,
#if 0
		"mov	UR,-(sp)\nmov	AR,-(sp)\n"
		"mov	UL,-(sp)\nmov	AL,-(sp)\n"
		"jsr	pc,ldiv\nadd	$10,sp\n",
#endif
		"ZK", },

{ DIV,	INCREG,
	SCREG,	TFLOAT,
	SCREG,	TFLOAT,
		0,	RLEFT,
		"divf	AR,AL\n", },

{ DIV,	INCREG,
	SCREG,			TDOUBLE,
	SCREG|SNAME|SOREG,	TDOUBLE,
		0,	RLEFT,
		"divf	AR,AL\n", },

/* XXX merge the two below to one */
{ MOD,	INBREG,
	SBREG|SNAME|SCON|SOREG,		TLONG,
	SBREG|SNAME|SCON|SOREG,		TLONG,
		NEEDS(NREG(B, 1), NSL(B), NSR(B), NEVER(R0),
		    NEVER(R1), NRES(R01)),	RESC1,
		"mov	UR,-(sp)\nmov	AR,-(sp)\n"
		"mov	UL,-(sp)\nmov	AL,-(sp)\n"
		"jsr	pc,lrem\nadd	$10,sp\n", },

{ MOD,	INBREG,
	SBREG|SNAME|SCON|SOREG,		TULONG,
	SBREG|SNAME|SCON|SOREG,		TULONG,
		NEEDS(NREG(B, 1), NSL(B), NSR(B), NEVER(R0),
		    NEVER(R1), NRES(R01)),	RESC1,
		"mov	UR,-(sp)\nmov	AR,-(sp)\n"
		"mov	UL,-(sp)\nmov	AL,-(sp)\n"
		"jsr	pc,ulrem\nadd	$10,sp\n", },

/* need extra move to be sure N flag is correct for sxt */
{ MOD,	INAREG,
	SAREG,		TINT|TPOINT,
	SAREG,		TINT|TPOINT,
		NEEDS(NEVER(R0), NEVER(R1), NLEFT(R1),
		    NRES(R1), NORIGHT(R0), NORIGHT(R1)), 	RDEST,
		"mov	AL,r1\nsxt	r0\ndiv	AR,r0\n", },

/* urem uses args in registers */
{ MOD,	INAREG,
	SAREG,		TUNSIGNED,
	SAREG,		TUNSIGNED,
		NEEDS(NREG(A, 1), NSL(A), NSR(A), NLEFT(R0), 
		    NRIGHT(R1), NRES(R0)),	RESC1,
		"jsr	pc,urem\n", },

/*
 * Indirection operators.
 */
{ UMUL,	INBREG,
	SANY,	TPOINT|TWORD,
	SOREG,	TLONG|TULONG,
		NBREG,	RESC1, /* |NBSL - may overwrite index reg */
		"mov	AL,A1\nmov	UL,U1\n", },

{ UMUL,	INAREG,
	SANY,	TPOINT|TWORD,
	SOREG,	TPOINT|TWORD,
		XSL(A),	RESC1,
		"mov	AL,A1\n", },

{ UMUL,	INAREG,
	SANY,	TANY,
	SOREG,	TCHAR|TUCHAR,
		XSL(A),	RESC1,
		"movb	AL,A1\n", },

{ UMUL,	INCREG,
	SANY,	TANY,
	SOREG,	TDOUBLE,
		NCREG,	RESC1,
		"movf	AL,A1\n", },

{ UMUL,	INCREG,
	SANY,	TANY,
	SOREG,	TFLOAT,
		NCREG,	RESC1,
		"movof	AL,A1\n", },

/*
 * Logical/branching operators
 */
{ OPLOG,	FORCC,
	SAREG|SOREG|SNAME|SCON,	TWORD|TPOINT,
	SZERO,	TANY,
		0, 	RESCC,
		"tst	AL\n", },

{ OPLOG,	FORCC,
	SAREG|SOREG|SNAME|SCON,	TCHAR|TUCHAR,
	SZERO,	TANY,
		0, 	RESCC,
		"tstb	AL\n", },

{ OPLOG,	FORCC,
	SAREG|SOREG|SNAME|SCON,	TWORD|TPOINT,
	SAREG|SOREG|SNAME|SCON,	TWORD|TPOINT,
		0, 	RESCC,
		"cmp	AL,AR\n", },

{ OPLOG,	FORCC,
	SCREG|SCON,	TFLOAT,
	SCREG,		TFLOAT,
		0, 	RESCC,
		"cmpf	AL,AR\ncfcc\n", },

{ OPLOG,	FORCC,
	SCREG|SOREG|SNAME|SCON,	TDOUBLE,
	SCREG,			TDOUBLE,
		0, 	RESCC,
		"cmpf	AL,AR\ncfcc\n", },

{ OPLOG,	FORCC,
	SAREG|SCON,	TCHAR|TUCHAR,
	SAREG|SCON,	TCHAR|TUCHAR,
		0, 	RESCC,
		"cmp	AL,AR\n", },

{ OPLOG,	FORCC,
	SBREG|SOREG|SNAME,	TLONG|TULONG,
	SBREG|SOREG|SNAME,	TLONG|TULONG,
		0,	RNULL,
		"ZF", },

{ AND,	INBREG|FORCC,
	SBREG,			TLONG|TULONG,
	SCON|SBREG|SOREG|SNAME,	TLONG|TULONG,
		0,	RLEFT|RESCC,
		"bic	AR,AL\nbic	UR,UL\n", },

/* set status bits */
{ AND,	FORCC,
	ARONS|SCON,	TWORD|TPOINT,
	ARONS|SCON,	TWORD|TPOINT,
		0,	RESCC,
		"bit	AR,AL\n", },

/* AND with int */
{ AND,	INAREG|FORCC|FOREFF,
	SAREG|SNAME|SOREG,	TWORD,
	SCON|SAREG|SOREG|SNAME,	TWORD,
		0,	RLEFT|RESCC,
		"bic	AR,AL\n", },

/* AND with char */
{ AND,	INAREG|FORCC,
	SAREG|SOREG|SNAME,	TCHAR|TUCHAR,
	ARONS|SCON,		TCHAR|TUCHAR,
		0,	RLEFT|RESCC,
		"bicb	AR,AL\n", },

{ OR,	INBREG|FORCC,
	SBREG,			TLONG|TULONG,
	SCON|SBREG|SOREG|SNAME,	TLONG|TULONG,
		0,	RLEFT|RESCC,
		"bis	AR,AL\nbis	UR,UL\n", },

/* OR with int */
{ OR,	FOREFF|INAREG|FORCC,
	ARONS,		TWORD,
	ARONS|SCON,	TWORD,
		0,	RLEFT|RESCC,
		"bis	AR,AL\n", },

/* OR with char */
{ OR,	INAREG|FORCC,
	SAREG|SOREG|SNAME,	TCHAR|TUCHAR,
	ARONS|SCON,		TCHAR|TUCHAR,
		0,	RLEFT|RESCC,
		"bisb	AR,AL\n", },

/* XOR with int (extended insn)  */
{ ER,	INAREG|FORCC,
	ARONS,	TWORD,
	SAREG,	TWORD,
		0,	RLEFT|RESCC,
		"xor	AR,AL\n", },

/* XOR with char (extended insn)  */
{ ER,	INBREG|FORCC,
	SBREG|SOREG|SNAME,	TLONG|TULONG,
	SBREG,	TLONG|TULONG,
		0,	RLEFT|RESCC,
		"xor	AR,AL\nxor	UR,UL\n", },

/*
 * Jumps.
 */
{ GOTO, 	FOREFF,
	SCON,	TANY,
	SANY,	TANY,
		0,	RNOP,
		"jbr	LL\n", },

{ GOTO, 	FOREFF,
	SAREG,	TANY,
	SANY,	TANY,
		0,	RNOP,
		"jmp	*AL\n", },

/*
 * Convert LTYPE to reg.
 */
/* Two bytes less if high half of constant is zero */
{ OPLTYPE,	INBREG,
	SANY,	TANY,
	SSCON,	TLONG|TULONG,
		NBREG,	RESC1,
		"mov	UL,U1\nsxt	A1\n", },

/* XXX - avoid OREG index register to be overwritten */
{ OPLTYPE,	INBREG,
	SANY,	TANY,
	SCON|SBREG|SNAME|SOREG,	TLONG|TULONG,
		NBREG,	RESC1,
		"mov	AL,A1\nmov	UL,U1\n", },

{ OPLTYPE,	INAREG,
	SANY,	TANY,
	SAREG|SCON|SOREG|SNAME,	TWORD|TPOINT,
		NEEDS(NREG(A, 1), NSR(A)), RESC1,
		"mov	AL,A1\n", },

{ OPLTYPE,	INAREG,
	SANY,	TANY,
	SAREG|SCON|SOREG|SNAME,	TCHAR,
		NAREG,		RESC1,
		"movb	AR,A1\n", },

{ OPLTYPE,	INAREG,
	SANY,	TANY,
	SAREG|SCON|SOREG|SNAME,	TUCHAR,
		NAREG,		RESC1,
		"clr	A1\nbisb	AL,A1\n", },

{ OPLTYPE,	INCREG,
	SANY,	TANY,
	SCREG|SCON|SOREG|SNAME,	TDOUBLE,
		NCREG,		RESC1,
		"movf	AL,A1\n", },

{ OPLTYPE,	INCREG,
	SANY,	TANY,
	SCREG|SCON|SOREG|SNAME,	TFLOAT,
		NCREG,		RESC1,
		"movof	AL,A1\n", },

/*
 * Negate a word.
 */
{ UMINUS,	INAREG|FOREFF,
	SAREG,	TWORD|TPOINT|TCHAR|TUCHAR,
	SANY,	TANY,
		0,	RLEFT,
		"neg	AL\n", },

{ UMINUS,	INBREG|FOREFF,
	SBREG|SOREG|SNAME,	TLONG,
	SANY,			TANY,
		0,	RLEFT,
		"neg	AL\nneg	UL\nsbc	AL\n", },

{ UMINUS,	INCREG|FOREFF,
	SCREG,	TFLOAT|TDOUBLE,
	SANY,	TANY,
		0,	RLEFT,
		"negf	AL\n", },


{ COMPL,	INBREG,
	SBREG,	TLONG|TULONG,
	SANY,	TANY,
		0,	RLEFT,
		"com	AL\ncom	UL\n", },

{ COMPL,	INAREG,
	SAREG,	TWORD,
	SANY,	TANY,
		0,	RLEFT,
		"com	AL\n", },

/*
 * Arguments to functions.
 */
{ FUNARG,	FOREFF,
	SBREG|SNAME|SOREG,	TLONG|TULONG,
	SANY,	TLONG|TULONG,
		0,	RNULL,
		"mov	UL,ZA(sp)\nmov	AL,-(sp)\n", },

{ FUNARG,	FOREFF,
	SCON,	TLONG|TULONG,
	SANY,	TLONG|TULONG,
		0,	RNULL,
		"ZD", },

{ FUNARG,	FOREFF,
	SZERO,	TANY,
	SANY,	TANY,
		0,	RNULL,
		"clr	ZA(sp)\n", },

{ FUNARG,	FOREFF,
	SARGSUB,	TWORD|TPOINT,
	SANY,		TWORD|TPOINT,
		0,	RNULL,
		"ZB", },

{ FUNARG,	FOREFF,
	SARGINC,	TWORD|TPOINT,
	SANY,		TWORD|TPOINT,
		0,	RNULL,
		"ZH", },

{ FUNARG,	FOREFF,
	SCON|SAREG|SNAME|SOREG,	TWORD|TPOINT,
	SANY,	TWORD|TPOINT,
		0,	RNULL,
		"mov	AL,ZA(sp)\n", },

{ FUNARG,	FOREFF,
	SCON,	TCHAR|TUCHAR,
	SANY,	TANY,
		0,	RNULL,
		"mov	AL,ZA(sp)\n", },

{ FUNARG,	FOREFF,
	SNAME|SOREG,	TCHAR,
	SANY,		TCHAR,
		NAREG,	RNULL,
		"movb	AL,A1\nmov	A1,ZA(sp)\n", },

{ FUNARG,	FOREFF,
	SNAME|SOREG,	TUCHAR,
	SANY,		TUCHAR,
		NAREG,	RNULL,
		"clr	ZA(sp)\nbisb	AL,(sp)\n", },

{ FUNARG,	FOREFF,
	SAREG,	TUCHAR|TCHAR,
	SANY,	TUCHAR|TCHAR,
		0,	RNULL,
		"mov	AL,ZA(sp)\n", },

{ FUNARG,	FOREFF,
	SCREG,	TFLOAT|TDOUBLE,
	SANY,		TANY,
		0,	RNULL,
		"movf	AL,ZA(sp)\n", },

{ STARG,	FOREFF,
	SAREG,	TPTRTO|TANY,
	SANY,	TSTRUCT,
		NEEDS(NEVER(R0), NLEFT(R1), NEVER(R2)),	0,
		"ZJ", },


# define DF(x) FORREW,SANY,TANY,SANY,TANY,NEEDS(NREWRITE),x,""

{ UMUL, DF( UMUL ), },

{ ASSIGN, DF(ASSIGN), },

{ STASG, DF(STASG), },

{ FLD, DF(FLD), },

{ OPLEAF, DF(NAME), },

/* { INIT, DF(INIT), }, */

{ OPUNARY, DF(UMINUS), },

{ OPANY, DF(BITYPE), },

{ FREE,	FREE,	FREE,	FREE,	FREE,	FREE,	0,	FREE,	"help; I'm in trouble\n" },
};

int tablesize = sizeof(table)/sizeof(table[0]);
