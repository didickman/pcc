#if 0
static char *sccsid ="@(#)table.c	1.33 (Berkeley) 5/11/88";
#endif

# include "pass2.h"

# define WPTR TPTRTO|TINT|TLONG|TFLOAT|TDOUBLE|TUNSIGNED|TULONG
# define TLL TLONGLONG|TULONGLONG
# define AWD SNAME|SOREG|SCON|STARNM|STARREG
/* tbl */
# define ANYSIGNED TINT|TLONG|TSHORT|TCHAR
# define ANYUSIGNED TUNSIGNED|TULONG|TUSHORT|TUCHAR
# define ANYFIXED ANYSIGNED|ANYUSIGNED
# define TUWORD TUNSIGNED|TULONG
# define TSWORD TINT|TLONG
# define TWORD TUWORD|TSWORD
# define NIAWD SNAME|SCON|STARNM
/* tbl */

struct optab table[] = {

#if 0
	/* the following entry is to fix a problem with
	   the manner that the first pass handles the
	   type of a shift expression                 */
{ PCONV,	INAREG|INTAREG,
	SAREG|AWD,	TINT|TUNSIGNED,
	SANY,	TPOINT,
		NAREG|NASL,	RLEFT,
		"", },

{ STASG,	INAREG,
	SNAME|SOREG,	TANY,
	SCON,	TANY,
		NAREG,	RESC1,
		"ZS	movl	AR,A1\n", },
#endif

/* Convert char pointer to int */
{ SCONV,	INTAREG,
	SAREG|STAREG,	TPTRTO|TCHAR|TUCHAR,
	SANY,	TWORD,
		NAREG,	RLEFT,
		"	lsh AL,2\n"
		"	move A1,AL\n"
		"	tlz AL,0740000\n"
		"	lsh A1,-042\n"
		"	trz A1,074\n"
		"	add AL,A1\n", },

/* Convert int to char pointer */
{ PCONV,	INTAREG,
	SAREG|STAREG,	TWORD,
	SANY,	TPTRTO|TCHAR|TUCHAR,
		NAREG,	RLEFT,
		"	move A1,AL\n"
		"	lsh A1,036\n"
		"	tlo A1,0700000\n"
		"	lsh AL,-2\n"
		"	ior AL,A1\n", },

/* Convert int/struct/foo pointer to char ptr */
{ PCONV,	INTAREG,
	STAREG,	TPOINT|TWORD|TSTRUCT,
	SANY,	TPTRTO|TCHAR|TUCHAR,
		0,	RLEFT,
		"	tlo AL,0700000\n", },

/* Convert int/struct/foo pointer to short ptr */
{ PCONV,	INTAREG,
	STAREG,	TPTRTO|TWORD|TSTRUCT,
	SANY,	TPTRTO|TSHORT|TUSHORT,
		0,	RLEFT,
		"	tlo AL,0740000\n", },

/* Convert char pointer to int/struct/multiple ptr */
{ PCONV,	INTAREG,
	STAREG,	TPTRTO|TCHAR|TUCHAR,
	SANY,	TPOINT|TWORD|TSTRUCT,
		0,	RLEFT,
		"	tlz AL,0770000\n", },

/* convert short/char to int. This is done when register is loaded */
{ SCONV,	INTAREG,
	STAREG,	TSHORT|TUSHORT|TCHAR|TUCHAR|TWORD,
	SANY,	TWORD,
		0,	RLEFT,
		"", },

/* convert int to short/char. This is done when register is loaded */
{ SCONV,	INTAREG,
	STAREG,	TWORD,
	SANY,	TSHORT|TUSHORT|TCHAR|TUCHAR|TWORD,
		0,	RLEFT,
		"", },

#if 0
{ SCONV,	INTAREG,
	SAREG|STAREG|SNAME|SOREG,	TSHORT,
	SANY,	TINT,
		NAREG|NASL,	RESC1,
		"	foo A1,AL\n", },
#endif

{ SCONV,	INTAREG,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SANY,	TULONGLONG,
		NAREG|NASL,	RESC1|RESC2,
		"	move A1,AL\n"
		"	setz U1\n", },

{ SCONV,	INTAREG,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SANY,	TLONGLONG,
		NAREG|NASL,	RESC1|RESC2,
		"	move A1,AL\n"
		"	move U1,A1\n"
		"	ash U1,-043\n", },

{ SCONV,	INTAREG,
	SAREG|STAREG|SNAME|SOREG,	TLL,
	SANY,	TWORD,
		NAREG|NASL,	RESC1,
		"	move A1,AL\n", },
#ifdef notyet
{ SCONV,	INTAREG|FORCC,
	SAREG|AWD,	TANY,
	SANY,	TANY,
		NAREG|NASL,	RESC1|RESCC,
		"	ZL\n", },

{ SCONV,	FORARG,
	SAREG|AWD,	TANY,
	SANY,	TANY,
		NAREG|NASL,	RNULL,
		"	ZM\n", },
#endif

#if 0
{ PMCONV,	INTAREG,
	STAREG,		TWORD,
	SCON,		TANY,
		NASL,		RESC1,
		"	adjbp AL,[ .long CR ]\n", },
#endif

/*
 * Store constant initializers.
 */
{ INIT,	FOREFF,
	SCON,	TANY,
	SANY,	TWORD|TPOINT,
		0,	RNOP,
		"	.long	CL\n", },

/*
 * Subroutine calls.
 */
{ UNARY CALL,	INTAREG,
	SCON,	TANY,
	SANY,	TWORD|TCHAR|TUCHAR|TSHORT|TUSHORT|TFLOAT|TDOUBLE|TLL|TPOINT,
		NAREG|NASL,     RESC1,
		"	pushj 017,ZI\n", },

{ UNARY CALL,	INTAREG,
	SAREG|STAREG|SNAME|SOREG,	TANY,
	SANY,	TWORD|TCHAR|TUCHAR|TSHORT|TUSHORT|TFLOAT|TDOUBLE|TLL|TPOINT,
		NAREG|NASL,	RESC1,	/* should be 0 */
		"	pushj 017,AL\n", },

/*
 * The next rules handle all "+="-style operators.
 */
{ ASG ER,	INAREG|FOREFF,
	STAREG|SAREG,	TWORD,
	SCON,		TWORD,
		0,	RLEFT,
		"ZS", },

{ ASG ER,	INAREG|FOREFF,
	STAREG|SAREG,	TLL,
	STAREG|SAREG,	TLL,
		0,	RLEFT,
		"	xor AL,AR\n"
		"	xor UL,UR\n", },

{ ASG OPSIMP,	INAREG|FOREFF,
	SAREG|STAREG,		TWORD,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
		0,	RLEFT,
		"	OR AL,AR\n", },

{ ASG OPSIMP,	INAREG|FOREFF,
	STAREG|SAREG,	TWORD,
	SCON,		TWORD,
		0,	RLEFT,
		"	ZF AL,ZG\n", },

{ ASG OPSIMP,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SAREG|STAREG,		TWORD,
		0,	RLEFT,
		"	OM AR,AL\n", },

{ ASG PLUS,	INAREG|FOREFF,
	SAREG|STAREG,			TLL,
	SAREG|STAREG|SNAME|SOREG,	TLL,
		0,	RLEFT,
		"	dadd AL,AR	# XXX - not accurate\n", },

/* Add to char/short pointer. XXX - should be able to remove the movem */
{ ASG PLUS,	INAREG|FOREFF,
	SAREG|STAREG|SOREG|SNAME,	TPTRTO|TCHAR|TUCHAR|TSHORT|TUSHORT,
	SAREG|STAREG,			TWORD,
		0,	RRIGHT,
		"	adjbp AR,AL\n"
		"	movem AR,AL\n", },

{ ASG PLUS,	INAREG|FOREFF,
	SAREG|STAREG,	TPTRTO|TCHAR|TUCHAR,
	SCON,		TWORD,
		0,	RLEFT,
		"ZX", },

{ ASG PLUS,     INAREG|FOREFF,
	SAREG|STAREG,			TWORD|TPOINT,
	SAREG|STAREG|SNAME|SOREG,	TWORD|TPOINT,
		0,	RLEFT,
		"	add AL,AR\n", },

{ ASG MINUS,     INAREG|FOREFF,
	SAREG|STAREG,			TWORD|TPOINT,
	SAREG|STAREG|SNAME|SOREG,	TWORD|TPOINT,
		0,	RLEFT,
		"	sub AL,AR\n", },

{ PLUS,	INAREG|FOREFF,
	SAREG|STAREG,	TPTRTO|TCHAR|TUCHAR,
	SCON,		TWORD,
		NAREG,	RESC1,
		"ZY", },

{ ASG OPSIMP,	INAREG|FOREFF,
	STAREG|SAREG,	TWORD|TPOINT,
	SCON,		TWORD,
		0,	RLEFT,
		"	ZF AL,ZG\n", },

{ ASG AND,	INAREG|FOREFF,
	SAREG|STAREG,			TLL,
	SAREG|STAREG|SNAME|SOREG,	TLL,
		0,	RLEFT,
		"	and AL,AR\n"
		"	and UL,UR\n", },


/*
 * The next rules handle all shift operators.
 */
{ ASG LS,	INTAREG|INAREG|FOREFF,
	SAREG|STAREG,	TWORD,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
		0,	RLEFT,
		"	OR AL,@AR\n", },

{ ASG LS,	INTAREG|INAREG|FOREFF,
	STAREG|SAREG,	TWORD,
	SCON,		TWORD,
		0,	RLEFT,
		"	ZF AL,ZH\n", },

{ ASG RS,	INTAREG|INAREG|FOREFF,
	STAREG|SAREG,	TSWORD,
	SCON,		TWORD,
		0,	RLEFT,
		"	ash AL,-ZH\n", },

{ ASG RS,	INTAREG|INAREG|FOREFF,
	STAREG|SAREG,	TUWORD,
	SCON,		TWORD,
		0,	RLEFT,
		"	lsh AL,-ZH\n", },

{ ASG LS,	INTAREG|INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SAREG|STAREG,		TWORD,
		0,	RLEFT,
		"	OM AR,@AL\n", },

{ ASG LS,       INTAREG|INAREG|FOREFF,
	STAREG|SAREG,	TLL,
	SCON,		TANY,
		0,	RLEFT,
		"	lshc AL,ZH\n", },

/*
 * The next rules takes care of assignments. "=".
 */
{ ASSIGN,	INAREG|INTAREG|FOREFF,
	SAREG,		TWORD,
	SAREG|SNAME|SOREG,	TWORD,
		0,	RLEFT,
		"	move AL,AR\n", },

{ ASSIGN,	INAREG|INTAREG|FOREFF,
	STAREG|SAREG,		TWORD,
	SCON,		TWORD,
		0,	RLEFT,
		"	movei AL,AR ZC\n", },

{ ASSIGN,	INAREG|INTAREG|FOREFF,
	SAREG|SNAME|SOREG,	TWORD,
	SAREG,		TWORD,
		0,	RRIGHT,
		"	movem AR,AL\n", },

{ ASSIGN,	INAREG|INTAREG|FOREFF,
	SAREG|SNAME|SOREG,	TWORD|TPOINT,
	SAREG|STAREG,		TWORD|TPOINT,
		0,	RRIGHT,
		"	movem AR,AL\n", },

{ ASSIGN,	INAREG|INTAREG|FOREFF,
	SAREG|STAREG,		TWORD,
	SAREG|SNAME|SOREG,	TWORD,
		0,	RLEFT,
		"	move AL,AR\n", },


{ ASSIGN,	INAREG|INTAREG|FOREFF,
	SAREG|SNAME|SOREG,	TLL,
	SAREG|STAREG,		TUWORD,
		0,	RLEFT,
		"	foomovem AR,AL\n", },

{ ASSIGN,	INAREG|INTAREG|FOREFF,
	SAREG|SNAME|SOREG,	TLL,
	SAREG|STAREG,		TLL,
		0,	RLEFT,
		"	dmovem AR,AL\n", },

{ ASSIGN,	INAREG|INTAREG|FOREFF,
	SOREG|SNAME,	TSHORT|TUSHORT|TCHAR|TUCHAR,
	SAREG|STAREG,	TANY,
		0,	RLEFT,
		"ZV", },

/*
 * DIV/MUL 
 * These can be done way more efficient.
 */
{ ASG DIV,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SCON,		TWORD,
		2*NAREG,	RLEFT,
		"	move A1,AL\n"
		"	setz U1,\n"
		"	idivi A1,AR\n"
		"	movem A1,AL\n", },

{ ASG DIV,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
		2*NAREG,	RLEFT,
		"	move A1,AL\n"
		"	setz U1,\n"
		"	idiv A1,AR\n"
		"	movem A1,AL\n", },

{ DIV,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SCON,		TWORD,
		2*NAREG,	RESC1,
		"	move A1,AL\n"
		"	setz U1,\n"
		"	idivi A1,AR\n", },

{ DIV,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
		2*NAREG,	RESC1,
		"	move A1,AL\n"
		"	setz U1,\n"
		"	idiv A1,AR\n", },

{ ASG MOD,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SCON,		TWORD,
		2*NAREG,	RLEFT,
		"	move A1,AL\n"
		"	setz U1,\n"
		"	idivi A1,AR\n"
		"	movem U1,AL\n", },

{ ASG MOD,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
		2*NAREG,	RLEFT,
		"	move A1,AL\n"
		"	setz U1,\n"
		"	idiv A1,AR\n"
		"	movem U1,AL\n", },

{ MOD,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SCON,		TWORD,
		2*NAREG,	RESC2,
		"	move A1,AL\n"
		"	setz U1,\n"
		"	idivi A1,AR\n", },

{ MOD,	INAREG|FOREFF,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
		2*NAREG,	RESC2,
		"	move A1,AL\n"
		"	setz U1,\n"
		"	idiv A1,AR\n", },

{ ASG MUL,	INAREG|FOREFF,
	SAREG|STAREG,	TWORD,
	SCON,		TWORD,
		0,	RLEFT,
		"Za", },

{ ASG MUL,	INAREG|FOREFF,
	SAREG|STAREG,			TWORD,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
		0,	RLEFT,
		"	imul AL,AR\n", },

/*
 * dummy UNARY MUL entry to get U* to possibly match OPLTYPE
 */
{ UNARY MUL,	FOREFF,
	SCC,	TANY,
	SCC,	TANY,
		0,	RNULL,
		"	HELP HELP HELP\n", },

/*
 * Logical/branching operators
 */
{ OPLTYPE,	FORCC,
	SANY,	TWORD|TPOINT,
	SANY,	TWORD|TPOINT,
		0,	RESCC,
		"	move 0,AR\n", },

{ OPLTYPE,	FORCC,
	SANY,	TLL,
	SANY,	TLL,
		0,	RESCC,
		"	move 0,AR\n	ior 0,UR\n", },

/* Match char/short pointers first, requires special handling */
{ OPLOG,	FORCC,
	SAREG|STAREG,	TPTRTO|TCHAR|TUCHAR|TSHORT|TUSHORT,
	SAREG|STAREG,	TPTRTO|TCHAR|TUCHAR|TSHORT|TUSHORT,
		0, 	RESCC,
		"ZZ", },

/* Can check anything by just comparing if EQ/NE */
{ EQ,		FORCC,
	SAREG|STAREG,	TWORD|TPOINT,
	SAREG|STAREG|SOREG|SNAME|SCON,	TWORD|TPOINT,
		0, 	RESCC,
		"ZR", },

{ NE,		FORCC,
	SAREG|STAREG,	TWORD|TPOINT,
	SAREG|STAREG|SOREG|SNAME|SCON,	TWORD|TPOINT,
		0, 	RESCC,
		"ZR", },

{ OPLOG,	FORCC,
	SAREG|STAREG,	TWORD,
	SAREG|STAREG|SOREG|SNAME|SCON,	TSWORD,
		0, 	RESCC,
		"ZR", },

{ OPLOG,	FORCC,  
	SAREG|STAREG,	TLL,
	SAREG|STAREG|SOREG|SNAME,	TLL,
		0,	RESCC,
		"ZQ", },

/*
 * Convert LTYPE to reg.
 */
{ OPLTYPE,	INAREG|INTAREG,
	SANY,	ANYFIXED,
	SCON,	ANYFIXED,
		NAREG|NASR,	RESC1,
		"	ZD A1,ZE\n", },

{ OPLTYPE,	INAREG|INTAREG,
	SANY,	TWORD|TPOINT,
	SCON,	TWORD|TPOINT,
		NAREG|NASR,	RESC1,
		"	xmovei A1,AR\n", },

{ OPLTYPE,	INAREG|INTAREG,
	SANY,	TWORD|TPOINT,
	SAREG|STAREG|SOREG|SNAME,	TWORD|TPOINT,
		NAREG|NASR,	RESC1,
		"	move A1,AR\n", },

{ OPLTYPE,	INAREG|INTAREG,
	SANY,	TLL,
	SCON,	TLL,
		NAREG,	RESC1,
		"	dmove A1,ZO\n", },

{ OPLTYPE,	INAREG|INTAREG,
	SANY,	TLL,
	SANY,	TLL,
		NAREG,	RESC1,
		"	dmove A1,AR\n", },

{ OPLTYPE,	INAREG|INTAREG,
	SOREG,		TSHORT|TUSHORT|TCHAR|TUCHAR,
	SANY,		TANY,
		NAREG|NASR,	RESC1,
		"ZU", },

/*
 * Negate a word.
 */
{ UNARY MINUS,	INAREG|INTAREG|FOREFF,
	SANY,	TWORD,
	SANY,	TWORD,
		NAREG|NASR,	RESC1,
		"	movn A1,AL\n", },

{ COMPL,	INTAREG,
	SAREG|STAREG|SNAME|SOREG,	TWORD,
	SANY,	TANY,
		NAREG|NASL,	RESC1|RESCC,
		"	setcam A1,AL\n", },

/*
 * Get condition codes.
 */
{ CCODES,	INAREG|INTAREG,
	SANY,	TANY,
	SANY,	TANY,
		NAREG,	RESC1,
		"	movei A1,01\nZN", },

/*
 * Arguments to functions.
 * These three should be possible to convert to one!
 */
{ REG,	FORARG,
	SANY,	TANY,
	SAREG|SNAME|SOREG,	TWORD|TPOINT,
		0,	RNULL,
		"	push 017,AR\n", },

{ OREG,	FORARG,
	SANY,	TANY,
	SAREG|SNAME|SOREG,	TWORD,
		0,	RNULL,
		"	push 017,AR\n", },

{ NAME,	FORARG,
	SANY,	TANY,
	SAREG|SNAME|SOREG,	TWORD,
		0,	RNULL,
		"	push 017,AR\n", },

{ ICON,	FORARG,
	SANY,	TANY,
	SCON,	TWORD,
		0,	RNULL,
		"	push 017,[ .long AR]\n", },

{ REG,	FORARG,
	SANY,		TANY,
	SAREG|STAREG,	TLL,
		0,	RNULL,
		"	push 017,AR\n	push 017,UR\n", },


# define DF(x) FORREW,SANY,TANY,SANY,TANY,REWRITE,x,""

{ UNARY MUL, DF( UNARY MUL ), },

{ INCR, DF(INCR), },

{ DECR, DF(INCR), },

{ ASSIGN, DF(ASSIGN), },

{ STASG, DF(STASG), },

{ FLD, DF(FLD), },

{ OPLEAF, DF(NAME), },

{ OPLOG,	FORCC,
	SANY,	TANY,
	SANY,	TANY,
		REWRITE,	BITYPE,
		"", },

{ OPLOG,	DF(NOT), },

{ COMOP, DF(COMOP), },

{ INIT, DF(INIT), },

{ OPUNARY, DF(UNARY MINUS), },

{ ASG OPANY, DF(ASG PLUS), },

{ OPANY, DF(BITYPE), },

{ FREE,	FREE,	FREE,	FREE,	FREE,	FREE,	FREE,	FREE,	"help; I'm in trouble\n" },
};
#if 0

struct optab  table[] = {

	/* the following entry is to fix a problem with
	   the manner that the first pass handles the
	   type of a shift expression                 */
{ PCONV,	INAREG|INTAREG,
	SAREG|AWD,	TINT|TUNSIGNED,
	SANY,	TPOINT,
		NAREG|NASL,	RLEFT,
		"", },

#if defined(FORT) || defined(SPRECC)
{ SCONV,	INTAREG|FORCC,
	SAREG|AWD,	TDOUBLE,
	SANY,	TFLOAT,
		NAREG|NASL,	RESC1|RESCC,
		"	cvtdf	AL,A1\n", },

{ SCONV,	INTAREG|FORCC,
	SAREG|AWD,	ANYSIGNED,
	SANY,	TFLOAT,
		NAREG|NASL,	RESC1|RESCC,
		"	cvtZLf	AL,TA1\n", },
#endif

/* take care of redundant conversions introduced by reclaim() */
{ SCONV,	INTAREG,
	STAREG,	TWORD,
	SANY,	TWORD,
		0,	RLEFT,
		"", },

{ SCONV,	INTAREG,
	STAREG,	TDOUBLE,
	SANY,	TDOUBLE,
		0,	RLEFT,
		"", },

{ SCONV,	INTAREG|FORCC,
	SAREG|AWD,	TANY,
	SANY,	TANY,
		NAREG|NASL,	RESC1|RESCC,
		"	ZA\n", },

{ SCONV,	FORARG,
	SAREG|AWD,	TANY,
	SANY,	TANY,
		NAREG|NASL,	RNULL,
		"	ZV\n", },

{ INIT,	FOREFF,
	SCON,	TANY,
	SANY,	TWORD,
		0,	RNOP,
		"	.long	CL\n", },

{ INIT,	FOREFF,
	SCON,	TANY,
	SANY,	TSHORT|TUSHORT,
		0,	RNOP,
		"	.word	CL\n", },

{ INIT,	FOREFF,
	SCON,	TANY,
	SANY,	TCHAR|TUCHAR,
		0,	RNOP,
		"	.byte	CL\n", },

#ifdef FORT
	/* for the use of fortran only */

{ GOTO,	FOREFF,
	SCON,	TANY,
	SANY,	TANY,
		0,	RNOP,
		"	jbr	CL\n", },
#endif

{ GOTO,	FOREFF,
	SNAME|SOREG,	TANY,
	SANY,	TANY,
		0,	RNOP,
		"	jmp	*AL\n", },

{ GOTO,	FOREFF,
	SAREG,	TANY,
	SANY,	TANY,
		0,	RNOP,
		"	jmp	(AL)\n", },

{ STARG,	FORARG,
	SCON|SOREG,	TANY,
	SANY,	TANY,
		0,	RNULL,
		"	subl2	ZT,sp\nZS", },

{ STASG,	FOREFF,
	SNAME|SOREG,	TANY,
	SCON|SAREG,	TANY,
		0,	RNOP,
		"ZS", },

{ STASG,	INAREG,
	SNAME|SOREG,	TANY,
	SCON,	TANY,
		NAREG,	RESC1,
		"ZS	movl	AR,A1\n", },

{ STASG,	INAREG,
	SNAME|SOREG,	TANY,
	SAREG,	TANY,
		0,	RRIGHT,
		"	pushl	AR\nZS	movl	(sp)+,AR\n", },

{ FLD,	INAREG|INTAREG,
	SANY,	TANY,
	SFLD,	ANYSIGNED,
		NAREG|NASR,	RESC1,
		"	extv	$H,$S,AR,A1\n", },

{ FLD,	INAREG|INTAREG,
	SANY,	TANY,
	SFLD,	ANYUSIGNED,
		NAREG|NASR,	RESC1,
		"	extzv	$H,$S,AR,A1\n", },

{ FLD,	FORARG,
	SANY,	TANY,
	SFLD,	ANYSIGNED,
		0,	RNULL,
		"	extv	$H,$S,AR,-(sp)\n", },

{ FLD,	FORARG,
	SANY,	TANY,
	SFLD,	ANYUSIGNED,
		0,	RNULL,
		"	extzv	$H,$S,AR,-(sp)\n", },

{ OPLOG,	FORCC,
	SAREG|AWD,	TWORD,
	SAREG|AWD,	TWORD,
		0,	RESCC,
		"	cmpl	AL,AR\nZP", },

{ OPLOG,	FORCC,
	SAREG|AWD,	TSHORT,
	SAREG|AWD,	TSHORT,
		0,	RESCC,
		"	cmpw	AL,AR\nZP", },

{ OPLOG,	FORCC,
	SAREG|AWD,	TUSHORT,
	SAREG|AWD,	TUSHORT,
		0,	RESCC,
		"	cmpw	AL,AR\nZP", },

{ OPLOG,	FORCC,
	SAREG|AWD,	TCHAR,
	SAREG|AWD,	TCHAR,
		0,	RESCC,
		"	cmpb	AL,AR\nZP", },

{ OPLOG,	FORCC,
	SAREG|AWD,	TUCHAR,
	SAREG|AWD,	TUCHAR,
		0,	RESCC,
		"	cmpb	AL,AR\nZP", },

/* optim2() handles degenerate comparisons with constants */
{ OPLOG,	FORCC,
	SAREG|AWD,	TCHAR|TUCHAR|TSHORT|TUSHORT,
	SCON,	ANYFIXED,
		0,	RESCC,
		"	cmpZL	AL,AR\nZP", },

{ OPLOG,	FORCC,
	SAREG|AWD,	TDOUBLE,
	SAREG|AWD,	TDOUBLE,
		0,	RESCC,
		"	cmpd	AL,AR\nZP", },

{ OPLOG,	FORCC,
	SAREG|AWD,	TFLOAT,
	SAREG|AWD,	TFLOAT,
		0,	RESCC,
		"	cmpf	AL,AR\nZP", },

#ifdef FORT
/* this really ought to be taken care of farther upstream... XXX */
{ OPLOG,	FORCC,
	SAREG|AWD,	TFLOAT,
	SAREG|AWD,	TDOUBLE,
		NAREG|NASL,	RESCC,
		"	cvtfd	AL,A1\n	cmpd	A1,AR\nZP", },

{ OPLOG,	FORCC,
	SAREG|AWD,	TDOUBLE,
	SAREG|AWD,	TFLOAT,
		NAREG|NASR,	RESCC,
		"	cvtfd	AR,A1\n	cmpd	AL,A1\nZP", },
#endif

{ CCODES,	INAREG|INTAREG,
	SANY,	TANY,
	SANY,	TANY,
		NAREG,	RESC1,
		"	movl	$1,A1\nZN", },

{ UNARY CALL,	INTAREG,
	SCON,	TANY,
	SANY,	TWORD|TCHAR|TUCHAR|TSHORT|TUSHORT|TFLOAT|TDOUBLE|TLL,
		NAREG|NASL,	RESC1,
		"	calls	ZC,CL\n", },

{ UNARY CALL,	INTAREG,
	SAREG,	TANY,
	SANY,	TWORD|TCHAR|TUCHAR|TSHORT|TUSHORT|TFLOAT|TDOUBLE|TLL,
		NAREG|NASL,	RESC1,	/* should be 0 */
		"	calls	ZC,(AL)\n", },

{ UNARY CALL,	INAREG|INTAREG,
	SNAME,	TANY,
	SANY,	TANY,
		NAREG|NASL,	RESC1,	/* really reg 0 */
		"	calls	ZC,*AL\n", },

{ UNARY CALL,	INAREG|INTAREG,
	SSOREG,	TANY,
	SANY,	TANY,
		NAREG|NASL,	RESC1,	/* really reg 0 */
		"	calls	ZC,*AL\n", },

{ ASG RS,	INAREG|FOREFF|FORCC,
	SAREG,	TWORD,
	SCON,	TINT|TUNSIGNED,
		0,	RLEFT|RESCC,
		"	extzv	AR,ZU,AL,AL\n", },

{ ASG RS,	INAREG|FOREFF|FORCC,
	SAREG,	TWORD,
	SAREG,	ANYFIXED,
		NAREG,	RLEFT|RESCC,
		"	subl3	AR,$32,A1\n	extzv	AR,A1,AL,AL\n", },

{ ASG RS,	INAREG|FOREFF|FORCC,
	SAREG,	TWORD,
	SAREG|AWD,	TWORD,
		NAREG,	RLEFT|RESCC,
		"	subl3	AR,$32,A1\n	extzv	AR,A1,AL,AL\n", },

{ RS,	INAREG|INTAREG|FORCC,
	SAREG,	TWORD,
	SCON,	TINT|TUNSIGNED,
		NAREG|NASL,	RESC1|RESCC,
		"	extzv	AR,ZU,AL,A1\n", },

{ ASG LS,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TWORD,
	SAREG|NIAWD,	ANYSIGNED|ANYUSIGNED,
		0,	RLEFT|RESCC,
		"	ashl	AR,AL,AL\n", },

{ ASG LS,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TWORD,
	SSOREG,	ANYSIGNED|ANYUSIGNED,
		0,	RLEFT|RESCC,
		"	ashl	AR,AL,AL\n", },

{ ASG LS,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TWORD,
	SOREG,	ANYSIGNED|ANYUSIGNED,
		NAREG,	RLEFT|RESCC,
		"	ZB	AR,A1\n	ashl	A1,AL,AL\n", },

/* long long addition */
{ ASG LS,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
	SAREG|NIAWD,	ANYSIGNED|ANYUSIGNED|TLONGLONG|TULONGLONG,
		0,	RLEFT|RESCC,
		"	ashq	AR,AL,AL\n", },

{ LS,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TWORD,
	SAREG|NIAWD,	ANYSIGNED|ANYUSIGNED,
		NAREG|NASL|NASR,	RESC1|RESCC,
		"	ashl	AR,AL,A1\n", },

{ LS,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TWORD,
	SSOREG,	ANYSIGNED|ANYUSIGNED,
		NAREG|NASL|NASR,	RESC1|RESCC,
		"	ashl	AR,AL,A1\n", },

{ LS,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TWORD,
	SOREG,	ANYSIGNED|ANYUSIGNED,
		NAREG|NASR,	RESC1|RESCC,
		"	ZB	AR,A1\n	ashl	A1,AL,A1\n", },

{ INCR,	FOREFF,
	SAREG|AWD,	TANY,
	SCON|SNAME,	TANY,
		0,	RLEFT,
		"	ZE\n", },

{ DECR,	FOREFF,
	SAREG|AWD,	TANY,
	SCON|SNAME,	TANY,
		0,	RLEFT,
		"	ZE\n", },

{ INCR,	INAREG|INTAREG,
	SAREG|AWD,	TANY,
	SCON|SNAME,	TANY,
		NAREG,	RESC1,
		"	ZD\n", },

{ DECR,	INAREG|INTAREG,
	SAREG|AWD,	TANY,
	SCON|SNAME,	TANY,
		NAREG,	RESC1,
		"	ZD\n", },

{ ASSIGN,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TFLOAT|TDOUBLE,
	SAREG|AWD,	TUCHAR|TUSHORT,
		NAREG|NASR,	RLEFT|RESCC,
		"	ZA\n", },

{ ASSIGN,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TANY,
	SAREG|AWD,	TANY,
		0,	RLEFT|RESCC,
		"	ZA\n", },

{ ASSIGN,	FOREFF,
	SFLD,	TANY,
	SAREG|AWD,	TWORD,
		0,	RNOP,
		"	insv	AR,$H,$S,AL\n", },

{ ASSIGN,	INAREG,
	SFLD,	ANYSIGNED,
	SAREG|AWD,	TWORD,
		NAREG,	RESC1,
		"	insv	AR,$H,$S,AL\n	extv	$H,$S,AL,A1\n", },

{ ASSIGN,	INAREG,
	SFLD,	ANYUSIGNED,
	SAREG|AWD,	TWORD,
		NAREG,	RESC1,
		"	insv	AR,$H,$S,AL\n	extzv	$H,$S,AL,A1\n", },

{ ASSIGN,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TWORD,
	SFLD,	ANYSIGNED,
		0,	RLEFT|RESCC,
		"	extv	$H,$S,AR,AL\n", },

{ ASSIGN,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TWORD,
	SFLD,	ANYUSIGNED,
		0,	RLEFT|RESCC,
		"	extzv	$H,$S,AR,AL\n", },

/* dummy UNARY MUL entry to get U* to possibly match OPLTYPE */
{ UNARY MUL,	FOREFF,
	SCC,	TANY,
	SCC,	TANY,
		0,	RNULL,
		"	HELP HELP HELP\n", },

{ OREG,	INTEMP,
	SANY,	TANY,
	SOREG,	TDOUBLE,
		2*NTEMP,	RESC1,
		"	movq	AR,A1\n", },

{ OREG,	INTEMP,
	SANY,	TANY,
	SOREG,	TANY,
		NTEMP,	RESC1,
		"	movZR	AR,A1\n", },

{ REG,	INTEMP,
	SANY,	TANY,
	SAREG,	TDOUBLE,
		2*NTEMP,	RESC1,
		"	movq	AR,A1\n", },

{ REG,	INTEMP,
	SANY,	TANY,
	SAREG,	TANY,
		NTEMP,	RESC1,
		"	movZF	AR,A1\n", },

#if defined(FORT) || defined(SPRECC)
{ REG,	FORARG,
	SANY,	TANY,
	SAREG,	TFLOAT,
		0,	RNULL,
		"	cvtfd	AR,-(sp)\n", },

{ REG,	FORARG,
	SANY,	TANY,
	SAREG,	TDOUBLE,
		0,	RNULL,
		"	movq	AR,-(sp)\n", },
#endif

{ OPLEAF,	FOREFF,
	SANY,	TANY,
	SAREG|AWD,	TANY,
		0,	RLEFT,
		"", },

{ OPLTYPE,	INAREG|INTAREG,
	SANY,	TANY,
	SANY,	TANY,
		NAREG|NASR,	RESC1,
		"	ZA\n", },

{ OPLTYPE,	FORCC,
	SANY,	TANY,
	SANY,	TANY,
		0,	RESCC,
		"	tstZR	AR\n", },

{ OPLTYPE,	FORARG,
	SANY,	TANY,
	SANY,	TANY,
		0,	RNULL,
		"	ZV\n", },

#if defined(FORT) || defined(SPRECC)
{ UNARY MINUS,	INTAREG|FORCC,
	SAREG|AWD,	TFLOAT,
	SANY,	TANY,
		NAREG|NASL,	RESC1|RESCC,
		"	mnegZL	TAL,A1\n", },

#endif

{ UNARY MINUS,	INTAREG|FORCC,
	SAREG|AWD,	TWORD|TDOUBLE,
	SANY,	TANY,
		NAREG|NASL,	RESC1|RESCC,
		"	mnegZL	AL,A1\n", },

{ COMPL,	INTAREG|FORCC,
	SAREG|AWD,	TWORD,
	SANY,	TANY,
		NAREG|NASL,	RESC1|RESCC,
		"	mcomZL	AL,A1\n", },

{ AND,	FORCC,
	SAREG|AWD,	TCHAR|TSHORT,
	SCON,	ANYFIXED,
		NAREG|NASL,	RESCC,
		"	ZZ\n", },

{ AND,	FORCC,
	SAREG|AWD,	TWORD|ANYUSIGNED,
	SCON,	ANYFIXED,
		0,	RESCC,
		"	ZZ\n", },

{ ASG AND,	INAREG|FOREFF|FORCC,
	SAREG,	TWORD,
	SCON,	TWORD,
		0,	RLEFT|RESCC,
		"	bicl2	AR,AL\n", },

/* General cases for DIV and ASG DIV are handled below with OPMUL */
/* Some special cases are handled in optim2() */

{ DIV,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TUNSIGNED|TULONG,
	SCON,	ANYUSIGNED,
		NAREG|NEVEN,	RESC1|RESCC,
		"	movl	AL,A1\n	clrl	U1\n	ediv	AR,A1,A1,U1\n", },

{ ASG DIV,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TINT|TLONG|TUNSIGNED|TULONG,
	SMCON,	ANYUSIGNED,
		0,	RLEFT|RESCC,
		"	ZJ\n", },

{ ASG DIV,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TINT|TLONG|TUNSIGNED|TULONG,
	SCON,	ANYUSIGNED,
		NAREG|NEVEN,	RLEFT|RESCC,
		"	movl	AL,A1\n	clrl	U1\n	ediv	AR,A1,AL,U1\n", },

{ MOD,	INAREG|INTAREG,
	SAREG|AWD,	TINT|TLONG,
	SAREG|AWD,	TINT|TLONG,
		NAREG,	RESC1,
		"	divl3	AR,AL,A1\n	mull2	AR,A1\n	subl3	A1,AL,A1\n", },

{ MOD,	INAREG|FOREFF,
	SAREG|AWD,	TUNSIGNED|TULONG,
	SMCON,	ANYUSIGNED,
		NAREG|NASL,	RLEFT|RESC1,
		"	ZJ\n", },

{ MOD,	INAREG|FOREFF,
	SAREG|AWD,	TUNSIGNED|TULONG,
	SCON,	ANYUSIGNED,
		NAREG|NEVEN,	RESC1,
		"	movl	AL,A1\n	clrl	U1\n	ediv	AR,A1,U1,A1\n", },

/* should only see UNSIGNED lhs here if converted from UCHAR/USHORT lhs */
{ ASG MOD,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TINT|TLONG|TUNSIGNED|TULONG,
	SAREG|AWD,	TINT|TLONG,
		NAREG,	RLEFT|RESCC,
		"	divl3	AR,AL,A1\n	mull2	AR,A1\n	subl2	A1,AL\n", },

{ ASG MOD,	INAREG|FOREFF,
	SAREG|AWD,	TINT|TLONG|TUNSIGNED|TULONG,
	SMCON,	ANYUSIGNED,
		0,	RLEFT,
		"	ZJ\n", },

{ ASG MOD,	INAREG|FOREFF,
	SAREG|AWD,	TINT|TLONG|TUNSIGNED|TULONG,
	SCON,	ANYUSIGNED,
		NAREG|NEVEN,	RLEFT,
		"	movl	AL,A1\n	clrl	U1\n	ediv	AR,A1,A1,AL\n", },

{ ASG OPMUL,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TINT|TUNSIGNED|TLONG|TULONG,
	SAREG|AWD,	TINT|TUNSIGNED|TLONG|TULONG,
		0,	RLEFT|RESCC,
		"	OL2	AR,AL\n", },

{ OPMUL,	INAREG|INTAREG|FORCC,
	STAREG,	TINT|TUNSIGNED|TLONG|TULONG,
	SAREG|AWD,	TINT|TUNSIGNED|TLONG|TULONG,
		0,	RLEFT|RESCC,
		"	OL2	AR,AL\n", },

{ OPMUL,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TINT|TUNSIGNED|TLONG|TULONG,
	SAREG|AWD,	TINT|TUNSIGNED|TLONG|TULONG,
		NAREG|NASL|NASR,	RESC1|RESCC,
		"	OL3	AR,AL,A1\n", },

{ ASG PLUS,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	ANYFIXED,
	SONE,	TANY,
		0,	RLEFT|RESCC,
		"	incZL	AL\n", },

{ ASG MINUS,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	ANYFIXED,
	SONE,	TANY,
		0,	RLEFT|RESCC,
		"	decZL	AL\n", },

{ PLUS,	INAREG|INTAREG|FORCC,
	STAREG,	TWORD,
	SONE,	TWORD,
		0,	RLEFT|RESCC,
		"	incZL	AL\n", },

{ MINUS,	INAREG|INTAREG|FORCC,
	STAREG,	TWORD,
	SONE,	TWORD,
		0,	RLEFT|RESCC,
		"	decZL	AL\n", },

{ ASG PLUS,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
		0,	RLEFT|RESCC,
		"	addl2\tAR,AL\n\tadwc\tUR,UL\n", },

{ ASG MINUS,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
		0,	RLEFT|RESCC,
		"	subl2\tAR,AL\n\tsbwc\tUR,UL\n", },

{ MUL,	INAREG|INTAREG|FOREFF,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
		2*NAREG,	RESC1|RESCC,
		"	emul	AL,AR,$0,A1	\n"
		"	ashl	$-31,AR,A2	\n"
		"	subl3	A2,UR,A2	\n"
		"	mull2	AL,A2		\n"
		"	addl2	A2,U1		\n"
		"	ashl	$-31,AL,A2	\n"
		"	subl3	A2,UL,A2	\n"
		"	mull2	AR,A2		\n"
		"	addl2	A2,U1		\n", },

/* XXX - rewriting rules should take care about this one */
{ ASG MUL,	INAREG|INTAREG|FOREFF,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
	SAREG|AWD,	TLONGLONG|TULONGLONG,
		2*NAREG,	RESC1|RESCC,
		"	emul	AL,AR,$0,A1	\n"
		"	ashl	$-31,AR,A2	\n"
		"	subl3	A2,UR,A2	\n"
		"	mull2	AL,A2		\n"
		"	addl2	A2,U1		\n"
		"	ashl	$-31,AL,A2	\n"
		"	subl3	A2,UL,A2	\n"
		"	mull2	AR,A2		\n"
		"	addl2	A2,U1		\n"
		"	movq	A1,AL		\n", },

/* Inefficient... */
{ ASG DIV,	INAREG|INTAREG|FOREFF,
	SAREG|AWD,	TLONGLONG|TULONGLONG, 
	SAREG|AWD,	TLONGLONG|TULONGLONG, 
		0,	RLEFT|RESCC,
		"	movq	AL,-(sp)\n"
		"	movq	AR,-(sp)\n"
		"	jsb	__pcc_divll\n"
		"	movq	(sp)+,AL\n", },

/* Inefficient... */
{ ASG MOD,	INAREG|INTAREG|FOREFF,
	SAREG|AWD,	TLONGLONG|TULONGLONG, 
	SAREG|AWD,	TLONGLONG|TULONGLONG, 
		0,	RLEFT|RESCC,
		"	movq	AL,-(sp)\n"
		"	movq	AR,-(sp)\n"
		"	jsb	__pcc_modll\n"
		"	movq	(sp)+,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TWORD,
	SAREG|AWD,	TWORD,
		0,	RLEFT|RESCC,
		"	OL2	AR,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	AWD,	TSHORT|TUSHORT,
	SAREG|SNAME|STARNM,	TSHORT|TUSHORT|TINT|TUNSIGNED|TLONG|TULONG,
		0,	RLEFT|RESCC,
		"	OW2	AR,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	AWD,	TSHORT|TUSHORT,
	SSOREG,	TSHORT|TUSHORT|TINT|TUNSIGNED|TLONG|TULONG,
		0,	RLEFT|RESCC,
		"	OW2	AR,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	AWD,	TSHORT|TUSHORT,
	SSCON,	ANYFIXED,
		0,	RLEFT|RESCC,
		"	OW2	AR,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	AWD,	TSHORT|TUSHORT,
	AWD,	TSHORT|TUSHORT,
		0,	RLEFT|RESCC,
		"	OW2	AR,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	AWD,	TCHAR|TUCHAR,
	SSOREG,	ANYFIXED,
		0,	RLEFT|RESCC,
		"	OB2	AR,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	AWD,	TCHAR|TUCHAR,
	SAREG|SNAME|STARNM,	ANYFIXED,
		0,	RLEFT|RESCC,
		"	OB2	AR,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	AWD,	TCHAR|TUCHAR,
	SCCON,	ANYFIXED,
		0,	RLEFT|RESCC,
		"	OB2	AR,AL\n", },

{ ASG OPSIMP,	INAREG|FOREFF|FORCC,
	AWD,	TCHAR|TUCHAR,
	AWD,	TCHAR|TUCHAR,
		0,	RLEFT|RESCC,
		"	OB2	AR,AL\n", },

{ OPSIMP,	INAREG|INTAREG|FORCC,
	STAREG,	ANYFIXED,
	SAREG|AWD,	TWORD,
		0,	RLEFT|RESCC,
		"	OL2	AR,AL\n", },

{ OPSIMP,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TWORD,
	SAREG|AWD,	TWORD,
		NAREG|NASL|NASR,	RESC1|RESCC,
		"	OL3	AR,AL,A1\n", },

{ ASG OPFLOAT,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TDOUBLE,
	SAREG|AWD,	TDOUBLE,
		0,	RLEFT|RESCC,
		"	OD2	AR,AL\n", },

{ ASG OPFLOAT,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	TFLOAT,
	SAREG|AWD,	TFLOAT,
		0,	RLEFT|RESCC,
#if defined(FORT) || defined(SPRECC)
		"	OF2	AR,TAL\n", },
#else
		"	OF2	AR,AL\n", },
#endif

{ ASG OPFLOAT,	INAREG|INTAREG|FOREFF|FORCC,
	SAREG|AWD,	TFLOAT,
	SAREG|AWD,	TDOUBLE,
		NAREG,	RLEFT|RESC1|RESCC,
		"	cvtfd	AL,A1\n	OD2	AR,A1\n	cvtdf	A1,AL\n", },

{ ASG OPFLOAT,	INAREG|FOREFF|FORCC,
	SAREG|AWD,	ANYFIXED,
#ifndef SPRECC
	SAREG|AWD,	TDOUBLE,		/* force FLOAT to register */
#else
	SAREG|AWD,	TFLOAT|TDOUBLE,
#endif
		NAREG,	RLEFT|RESCC,	/* usable() knows we need a reg pair */
		"	ZG\n", },

{ OPFLOAT,	INAREG|INTAREG|FORCC,
	STAREG,	TDOUBLE,
	SAREG|AWD,	TDOUBLE,
		0,	RLEFT|RESCC,
		"	OD2	AR,AL\n", },

{ OPFLOAT,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TDOUBLE,
	SAREG|AWD,	TDOUBLE,
		NAREG|NASL|NASR,	RESC1|RESCC,
		"	OD3	AR,AL,A1\n", },

#if defined(FORT) || defined(SPRECC)
{ OPFLOAT,	INAREG|INTAREG|FORCC,
	STAREG,		TFLOAT,
	SAREG|AWD,	TFLOAT,
		0,	RLEFT|RESCC,
		"	OF2	AR,TAL\n", },

{ OPFLOAT,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TFLOAT,
	SAREG|AWD,	TFLOAT,
		NAREG|NASL|NASR,	RESC1|RESCC,
		"	OF3	AR,AL,TA1\n", },
#endif

#ifdef FORT
/* perform some implicit conversions XXX SHOULD FIX f77 FRONT END */
{ OPFLOAT,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TFLOAT,
	SAREG|AWD,	TDOUBLE,
		NAREG|NASL,	RESC1|RESCC,
		"	cvtfd	AL,A1\n	OD2	AR,A1\n", },

{ OPFLOAT,	INAREG|INTAREG|FORCC,
	SAREG|AWD,	TDOUBLE,
	SAREG|AWD,	TFLOAT,
		NAREG|NASR,	RESC1|RESCC,
		"	cvtfd	AR,A1\n	OD3	A1,AL,A1\n", },
#endif

	/* Default actions for hard trees ... */

# define DF(x) FORREW,SANY,TANY,SANY,TANY,REWRITE,x,""

{ UNARY MUL, DF( UNARY MUL ), },

{ INCR, DF(INCR), },

{ DECR, DF(INCR), },

{ ASSIGN, DF(ASSIGN), },

{ STASG, DF(STASG), },

{ FLD, DF(FLD), },

{ OPLEAF, DF(NAME), },

{ OPLOG,	FORCC,
	SANY,	TANY,
	SANY,	TANY,
		REWRITE,	BITYPE,
		"", },

{ OPLOG,	DF(NOT), },

{ COMOP, DF(COMOP), },

{ INIT, DF(INIT), },

{ OPUNARY, DF(UNARY MINUS), },

{ ASG OPANY, DF(ASG PLUS), },

{ OPANY, DF(BITYPE), },

{ FREE,	FREE,	FREE,	FREE,	FREE,	FREE,	FREE,	FREE,	"help; I'm in trouble\n" },
};
#endif
