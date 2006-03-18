
struct optab;
struct symtab;
struct sw;

void cerror(char *s, ...);
void werror(char *s, ...);
void uerror(char *s, ...);
void reclaim(NODE *p, int, int);
void walkf(NODE *, void (*f)(NODE *));
void allchk(void);
void tfree(NODE *);
int tshape(NODE *, int);
void prtdcon(NODE *p);
void tinit(void);
void tcheck(void);
void mkdope(void);
int tshape(NODE *p, int shape);
int shtemp(NODE *p);
int flshape(NODE *p);
int shumul(NODE *p);
int ttype(TWORD t, int tword);
void expand(NODE *, int, char *);
void hopcode(int, int);
void adrcon(CONSZ);
void zzzcode(NODE *, int);
void insput(NODE *);
void upput(NODE *, int);
void econvert(NODE *);
int andable(NODE *);
int conval(NODE *, int, NODE *);
int ispow2(CONSZ);
void defid(NODE *q, int class);
int getlab(void);
void ftnend(void);
void efcode(void);
void dclargs(void);
void fixarg(struct symtab *);
void cendarg(void);
void defalign(int);
int fldal(unsigned int);
void vfdzero(int);
void zecode(int);
void ilbrace(void);
void irbrace(void);
void irbrace(void);
void putbyte(int v);
void endinit(void);
void doinit(NODE *p);
void ecomp(NODE *p);
void cinit(NODE *, int);
void bccode(void);
int upoff(int size, int alignment, int *poff);
void fldty(struct symtab *p);
void nidcl(NODE *p, int class);
int noinit(void);
void eprint(NODE *, int, int *, int *);
int uclass(int class);
int fixclass(int, TWORD type);
void lineid(int, char *);
void mycanon(NODE *);
void delay(NODE *);
int delay1(NODE *);
void delay2(NODE *);
void setregs(void);
int autoincr(NODE *);
int deltest(NODE *);
void canon(NODE *);
void order(NODE *, int);
int tlen(NODE *p);
int setincr(NODE *);
int setbin(NODE *);
void stoarg(NODE *p, int);
void constore(NODE *);
void markcall(NODE *);
void oreg2(NODE *p);
int notoff(TWORD, int, CONSZ, char *);
void bycode(int, int);
void pstab(char *, int);
void psline(void);
int notlval(NODE *);
int icons(NODE *);
void ecode(NODE *p);
int yylex(void);
void yyerror(char *s);
void p2tree(NODE *p);
int rewfld(NODE *p);
int freetemp(int k);
