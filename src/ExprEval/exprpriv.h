/*
    File: exprpriv.h
    Auth: Brian Allen Vanderburg II
    Date: Tuesday, February 28, 2006
    Desc: Private include file for ExprEval library

    This file is part of ExprEval.
*/


/* Include once */
#ifndef __BAVII_EXPRPRIV_H
#define __BAVII_EXPRPRIV_H

/* Need some definitions, NULL, etc */
#include <stddef.h>

/* Include config and main expreval header */
#include "expreval.h"
#include "exprconf.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Version number
*/
#define EXPR_VERSIONMAJOR 2
#define EXPR_VERSIONMINOR 7

/* Node types */
enum
    {
    EXPR_NODETYPE_UNKNOWN = 0,
    EXPR_NODETYPE_MULTI,
    EXPR_NODETYPE_ADD,
    EXPR_NODETYPE_SUBTRACT,
    EXPR_NODETYPE_MULTIPLY,
    EXPR_NODETYPE_DIVIDE,
    EXPR_NODETYPE_EXPONENT,
    EXPR_NODETYPE_NEGATE,
    EXPR_NODETYPE_VALUE,
    EXPR_NODETYPE_VARIABLE,
    EXPR_NODETYPE_ASSIGN,
    EXPR_NODETYPE_FUNCTION
    };

/* Functions can be evaluated directly in EXPREVAL.  If fptr
   is NULL, type is used to determine what the function is */
enum
    {
    EXPR_NODEFUNC_UNKNOWN = 0,
    EXPR_NODEFUNC_EXIT,
    EXPR_NODEFUNC_MSS,
    EXPR_NODEFUNC_ABS,
    EXPR_NODEFUNC_MOD,
    EXPR_NODEFUNC_IPART,
    EXPR_NODEFUNC_FPART,
    EXPR_NODEFUNC_MIN,
    EXPR_NODEFUNC_MAX,
    //
    EXPR_NODEFUNC_BITCOUNTER,

    EXPR_NODEFUNC_VERSION,
    EXPR_NODEFUNC_ALGEBRA,
    //
    EXPR_NODEFUNC_BINSUM,
    EXPR_NODEFUNC_BINSUB,
    EXPR_NODEFUNC_COMP,
    //
    EXPR_NODEFUNC_POW,
    EXPR_NODEFUNC_SQRT,
    //
    EXPR_NODEFUNC_CBRT,
    EXPR_NODEFUNC_ROOT,
    //
    EXPR_NODEFUNC_SIN,
    EXPR_NODEFUNC_SINH,
    //
    EXPR_NODEFUNC_CSC,
    EXPR_NODEFUNC_CSCH,
    EXPR_NODEFUNC_ASIN,
    EXPR_NODEFUNC_ASINH,
    EXPR_NODEFUNC_ACSC,
    EXPR_NODEFUNC_ACSCH,
    EXPR_NODEFUNC_COS,
    EXPR_NODEFUNC_COSH,
    EXPR_NODEFUNC_SEC,
    EXPR_NODEFUNC_SECH,
    EXPR_NODEFUNC_ACOS,
    EXPR_NODEFUNC_ACOSH,
    EXPR_NODEFUNC_ASEC,
    EXPR_NODEFUNC_ASECH,
    EXPR_NODEFUNC_TAN,
    EXPR_NODEFUNC_TANH,
    EXPR_NODEFUNC_COT,
    EXPR_NODEFUNC_COTH,
    EXPR_NODEFUNC_ATAN,
    EXPR_NODEFUNC_ATANH,
    EXPR_NODEFUNC_ACOT,
    EXPR_NODEFUNC_ACOTH,
    EXPR_NODEFUNC_HSIN,
    EXPR_NODEFUNC_HSINH,
    EXPR_NODEFUNC_QSIN,
    EXPR_NODEFUNC_QSINH,
    EXPR_NODEFUNC_HCOS,
    EXPR_NODEFUNC_HCOSH,
    EXPR_NODEFUNC_QCOS,
    EXPR_NODEFUNC_QCOSH,
    EXPR_NODEFUNC_HSEC,
    EXPR_NODEFUNC_HSECH,
    EXPR_NODEFUNC_QSEC,
    EXPR_NODEFUNC_QSECH,
    EXPR_NODEFUNC_HCSC,
    EXPR_NODEFUNC_HCSCH,
    EXPR_NODEFUNC_QCSC,
    EXPR_NODEFUNC_QCSCH,
    EXPR_NODEFUNC_HTAN,
    EXPR_NODEFUNC_HTANH,
    EXPR_NODEFUNC_QTAN,
    EXPR_NODEFUNC_QTANH,
    EXPR_NODEFUNC_HCOT,
    EXPR_NODEFUNC_HCOTH,
    EXPR_NODEFUNC_QCOT,
    EXPR_NODEFUNC_QCOTH,
    EXPR_NODEFUNC_VSIN,
    EXPR_NODEFUNC_VSINH,
    EXPR_NODEFUNC_CVSIN,
    EXPR_NODEFUNC_CVSINH,
    EXPR_NODEFUNC_VCOS,
    EXPR_NODEFUNC_VCOSH,
    EXPR_NODEFUNC_CVCOS,
    EXPR_NODEFUNC_CVCOSH,
    EXPR_NODEFUNC_HVSIN,
    EXPR_NODEFUNC_HVSINH,
    EXPR_NODEFUNC_HCVSIN,
    EXPR_NODEFUNC_HCVSINH,
    EXPR_NODEFUNC_QVSIN,
    EXPR_NODEFUNC_QVSINH,
    EXPR_NODEFUNC_QCVSIN,
    EXPR_NODEFUNC_QCVSINH,
    EXPR_NODEFUNC_HVCOS,
    EXPR_NODEFUNC_HVCOSH,
    EXPR_NODEFUNC_HCVCOS,
    EXPR_NODEFUNC_HCVCOSH,
    EXPR_NODEFUNC_QVCOS,
    EXPR_NODEFUNC_QVCOSH,
    EXPR_NODEFUNC_QCVCOS,
    EXPR_NODEFUNC_QCVCOSH,
    EXPR_NODEFUNC_ESEC,
    EXPR_NODEFUNC_ESECH,
    EXPR_NODEFUNC_ECSC,
    EXPR_NODEFUNC_ECSCH,
    EXPR_NODEFUNC_HESEC,
    EXPR_NODEFUNC_HESECH,
    EXPR_NODEFUNC_HECSC,
    EXPR_NODEFUNC_HECSCH,
    EXPR_NODEFUNC_QESEC,
    EXPR_NODEFUNC_QESECH,
    EXPR_NODEFUNC_QECSC,
    EXPR_NODEFUNC_QECSCH,
    EXPR_NODEFUNC_SINC,
    EXPR_NODEFUNC_SINCH,
    EXPR_NODEFUNC_HSINC,
    EXPR_NODEFUNC_HSINCH,
    EXPR_NODEFUNC_QSINC,
    EXPR_NODEFUNC_QSINCH,
    EXPR_NODEFUNC_COSC,
    EXPR_NODEFUNC_COSCH,
    EXPR_NODEFUNC_HCOSC,
    EXPR_NODEFUNC_HCOSCH,
    EXPR_NODEFUNC_QCOSC,
    EXPR_NODEFUNC_QCOSCH,
    EXPR_NODEFUNC_SECC,
    EXPR_NODEFUNC_SECCH,
    EXPR_NODEFUNC_HSECC,
    EXPR_NODEFUNC_HSECCH,
    EXPR_NODEFUNC_QSECC,
    EXPR_NODEFUNC_QSECCH,
    EXPR_NODEFUNC_CSCC,
    EXPR_NODEFUNC_CSCCH,
    EXPR_NODEFUNC_HCSCC,
    EXPR_NODEFUNC_HCSCCH,
    EXPR_NODEFUNC_QCSCC,
    EXPR_NODEFUNC_QCSCCH,
    EXPR_NODEFUNC_TANC,
    EXPR_NODEFUNC_TANCH,
    EXPR_NODEFUNC_HTANC,
    EXPR_NODEFUNC_HTANCH,
    EXPR_NODEFUNC_QTANC,
    EXPR_NODEFUNC_QTANCH,
    EXPR_NODEFUNC_COTC,
    EXPR_NODEFUNC_COTCH,
    EXPR_NODEFUNC_HCOTC,
    EXPR_NODEFUNC_HCOTCH,
    EXPR_NODEFUNC_QCOTC,
    EXPR_NODEFUNC_QCOTCH,
    //
    EXPR_NODEFUNC_ATAN2,
    
    EXPR_NODEFUNC_CSIN,
    EXPR_NODEFUNC_CSINH,
    //
    EXPR_NODEFUNC_CCSC,
    EXPR_NODEFUNC_CCSCH,
    EXPR_NODEFUNC_CASIN,
    EXPR_NODEFUNC_CASINH,
    EXPR_NODEFUNC_CACSC,
    EXPR_NODEFUNC_CACSCH,
    EXPR_NODEFUNC_CCOS,
    EXPR_NODEFUNC_CCOSH,
    EXPR_NODEFUNC_CSEC,
    EXPR_NODEFUNC_CSECH,
    EXPR_NODEFUNC_CACOS,
    EXPR_NODEFUNC_CACOSH,
    EXPR_NODEFUNC_CASEC,
    EXPR_NODEFUNC_CASECH,
    EXPR_NODEFUNC_CTAN,
    EXPR_NODEFUNC_CTANH,
    EXPR_NODEFUNC_CCOT,
    EXPR_NODEFUNC_CCOTH,
    EXPR_NODEFUNC_CATAN,
    EXPR_NODEFUNC_CATANH,
    EXPR_NODEFUNC_CACOT,
    EXPR_NODEFUNC_CACOTH,
    EXPR_NODEFUNC_CHSIN,
    EXPR_NODEFUNC_CHSINH,
    EXPR_NODEFUNC_CQSIN,
    EXPR_NODEFUNC_CQSINH,
    EXPR_NODEFUNC_CHCOS,
    EXPR_NODEFUNC_CHCOSH,
    EXPR_NODEFUNC_CQCOS,
    EXPR_NODEFUNC_CQCOSH,
    EXPR_NODEFUNC_CHSEC,
    EXPR_NODEFUNC_CHSECH,
    EXPR_NODEFUNC_CQSEC,
    EXPR_NODEFUNC_CQSECH,
    EXPR_NODEFUNC_CHCSC,
    EXPR_NODEFUNC_CHCSCH,
    EXPR_NODEFUNC_CQCSC,
    EXPR_NODEFUNC_CQCSCH,
    EXPR_NODEFUNC_CHTAN,
    EXPR_NODEFUNC_CHTANH,
    EXPR_NODEFUNC_CQTAN,
    EXPR_NODEFUNC_CQTANH,
    EXPR_NODEFUNC_CHCOT,
    EXPR_NODEFUNC_CHCOTH,
    EXPR_NODEFUNC_CQCOT,
    EXPR_NODEFUNC_CQCOTH,
    EXPR_NODEFUNC_CPXVSIN,
    EXPR_NODEFUNC_CPXVSINH,
    EXPR_NODEFUNC_CCVSIN,
    EXPR_NODEFUNC_CCVSINH,
    EXPR_NODEFUNC_CPXVCOS,
    EXPR_NODEFUNC_CPXVCOSH,
    EXPR_NODEFUNC_CCVCOS,
    EXPR_NODEFUNC_CCVCOSH,
    EXPR_NODEFUNC_CHVSIN,
    EXPR_NODEFUNC_CHVSINH,
    EXPR_NODEFUNC_CHCVSIN,
    EXPR_NODEFUNC_CHCVSINH,
    EXPR_NODEFUNC_CQVSIN,
    EXPR_NODEFUNC_CQVSINH,
    EXPR_NODEFUNC_CQCVSIN,
    EXPR_NODEFUNC_CQCVSINH,
    EXPR_NODEFUNC_CHVCOS,
    EXPR_NODEFUNC_CHVCOSH,
    EXPR_NODEFUNC_CHCVCOS,
    EXPR_NODEFUNC_CHCVCOSH,
    EXPR_NODEFUNC_CQVCOS,
    EXPR_NODEFUNC_CQVCOSH,
    EXPR_NODEFUNC_CQCVCOS,
    EXPR_NODEFUNC_CQCVCOSH,
    EXPR_NODEFUNC_CESEC,
    EXPR_NODEFUNC_CESECH,
    EXPR_NODEFUNC_CECSC,
    EXPR_NODEFUNC_CECSCH,
    EXPR_NODEFUNC_CHESEC,
    EXPR_NODEFUNC_CHESECH,
    EXPR_NODEFUNC_CHECSC,
    EXPR_NODEFUNC_CHECSCH,
    EXPR_NODEFUNC_CQESEC,
    EXPR_NODEFUNC_CQESECH,
    EXPR_NODEFUNC_CQECSC,
    EXPR_NODEFUNC_CQECSCH,
    EXPR_NODEFUNC_CSINC,
    EXPR_NODEFUNC_CSINCH,
    EXPR_NODEFUNC_CHSINC,
    EXPR_NODEFUNC_CHSINCH,
    EXPR_NODEFUNC_CQSINC,
    EXPR_NODEFUNC_CQSINCH,
    EXPR_NODEFUNC_CCOSC,
    EXPR_NODEFUNC_CCOSCH,
    EXPR_NODEFUNC_CHCOSC,
    EXPR_NODEFUNC_CHCOSCH,
    EXPR_NODEFUNC_CQCOSC,
    EXPR_NODEFUNC_CQCOSCH,
    EXPR_NODEFUNC_CSECC,
    EXPR_NODEFUNC_CSECCH,
    EXPR_NODEFUNC_CHSECC,
    EXPR_NODEFUNC_CHSECCH,
    EXPR_NODEFUNC_CQSECC,
    EXPR_NODEFUNC_CQSECCH,
    EXPR_NODEFUNC_CCSCC,
    EXPR_NODEFUNC_CCSCCH,
    EXPR_NODEFUNC_CHCSCC,
    EXPR_NODEFUNC_CHCSCCH,
    EXPR_NODEFUNC_CQCSCC,
    EXPR_NODEFUNC_CQCSCCH,
    EXPR_NODEFUNC_CTANC,
    EXPR_NODEFUNC_CTANCH,
    EXPR_NODEFUNC_CHTANC,
    EXPR_NODEFUNC_CHTANCH,
    EXPR_NODEFUNC_CQTANC,
    EXPR_NODEFUNC_CQTANCH,
    EXPR_NODEFUNC_CCOTC,
    EXPR_NODEFUNC_CCOTCH,
    EXPR_NODEFUNC_CHCOTC,
    EXPR_NODEFUNC_CHCOTCH,
    EXPR_NODEFUNC_CQCOTC,
    EXPR_NODEFUNC_CQCOTCH,
    //
    EXPR_NODEFUNC_MATRIXDET,
    EXPR_NODEFUNC_MATRIXNORM,
    EXPR_NODEFUNC_MATRIXTRACE,
    EXPR_NODEFUNC_MATRIXRANK,
    EXPR_NODEFUNC_MATRIXILLCHK,
    EXPR_NODEFUNC_SCALARPROD,
    //
    EXPR_NODEFUNC_LOG,
    EXPR_NODEFUNC_LOG2,
    EXPR_NODEFUNC_POW10,
    EXPR_NODEFUNC_LN,
    EXPR_NODEFUNC_EXP,
    EXPR_NODEFUNC_EXPC,
    EXPR_NODEFUNC_EXP10,
    EXPR_NODEFUNC_EXP10C,
    EXPR_NODEFUNC_EXP2,
    EXPR_NODEFUNC_EXP2C,
    EXPR_NODEFUNC_LOGN,
    EXPR_NODEFUNC_LOGC,
    EXPR_NODEFUNC_LNC,
    EXPR_NODEFUNC_LOG2C,
    EXPR_NODEFUNC_LOG1P,
    EXPR_NODEFUNC_LOG1PC,
    EXPR_NODEFUNC_LOG101P,
    EXPR_NODEFUNC_LOG101PC,
    EXPR_NODEFUNC_LOG21P,
    EXPR_NODEFUNC_LOG21PC,
    EXPR_NODEFUNC_CEXP,
    EXPR_NODEFUNC_CEXPC,
    EXPR_NODEFUNC_CEXP10,
    EXPR_NODEFUNC_CEXP10C,
    EXPR_NODEFUNC_CEXP2,
    EXPR_NODEFUNC_CEXP2C,
    EXPR_NODEFUNC_CPOW,
    EXPR_NODEFUNC_CROOT,
    EXPR_NODEFUNC_CSQRT,
    EXPR_NODEFUNC_CCBRT,
    EXPR_NODEFUNC_CLOGN,
    EXPR_NODEFUNC_CLN,
    EXPR_NODEFUNC_CLNC,
    EXPR_NODEFUNC_CLOG,
    EXPR_NODEFUNC_CLOGC,
    EXPR_NODEFUNC_CLOG2,
    EXPR_NODEFUNC_CLOG2C,
    EXPR_NODEFUNC_CLOG1P,
    EXPR_NODEFUNC_CLOG1PC,
    EXPR_NODEFUNC_CLOG101P,
    EXPR_NODEFUNC_CLOG101PC,
    EXPR_NODEFUNC_CLOG21P,
    EXPR_NODEFUNC_CLOG21PC,
    EXPR_NODEFUNC_CARG,
    EXPR_NODEFUNC_CABS,
    EXPR_NODEFUNC_QABS,
    EXPR_NODEFUNC_OABS,
    EXPR_NODEFUNC_SABS,
    EXPR_NODEFUNC_CEIL,
    EXPR_NODEFUNC_FLOOR,
    EXPR_NODEFUNC_SGEQSOLVER,
    EXPR_NODEFUNC_COMPLEXADD,
    EXPR_NODEFUNC_COMPLEXMUL,
    EXPR_NODEFUNC_QUATERNIONSADD,
    EXPR_NODEFUNC_QUATERNIONSMUL,
    EXPR_NODEFUNC_RAND,
    EXPR_NODEFUNC_RANDOM,
    EXPR_NODEFUNC_RANDOMIZE,
    EXPR_NODEFUNC_DEG,
    EXPR_NODEFUNC_RAD,
    EXPR_NODEFUNC_RECTTOPOLR,
    EXPR_NODEFUNC_RECTTOPOLA,
    EXPR_NODEFUNC_POLTORECTX,
    EXPR_NODEFUNC_POLTORECTY,
    EXPR_NODEFUNC_CBASE,
    EXPR_NODEFUNC_NPNUM,
    EXPR_NODEFUNC_PRIMORIAL,
    EXPR_NODEFUNC_FPNSUM,
    EXPR_NODEFUNC_FIBONACCIAL,
    EXPR_NODEFUNC_LCM,
    EXPR_NODEFUNC_GCD,
    EXPR_NODEFUNC_FACT,
    EXPR_NODEFUNC_SFACT,
    EXPR_NODEFUNC_STIRLING,
    EXPR_NODEFUNC_FIBO,
    EXPR_NODEFUNC_PERMS,
    EXPR_NODEFUNC_PERMSREP,
    EXPR_NODEFUNC_KPERMS,
    EXPR_NODEFUNC_KPERMSREP,
    EXPR_NODEFUNC_COMBS,
    EXPR_NODEFUNC_COMBSREP,
    EXPR_NODEFUNC_COMB,
    EXPR_NODEFUNC_GSUM,
    EXPR_NODEFUNC_ASUM,
    EXPR_NODEFUNC_GASUM,
    EXPR_NODEFUNC_FSUM,
    EXPR_NODEFUNC_FASUM,
    EXPR_NODEFUNC_SFASUM,
    EXPR_NODEFUNC_FNNSUM,
    //
    EXPR_NODEFUNC_SUM,
    EXPR_NODEFUNC_PRODUCT,
    EXPR_NODEFUNC_MEDIA,
    EXPR_NODEFUNC_MODE,
    EXPR_NODEFUNC_VARIANCE,
    EXPR_NODEFUNC_COVARIANCE,
    EXPR_NODEFUNC_STDDEV,
    EXPR_NODEFUNC_OUTLIER,
    EXPR_NODEFUNC_OUTLIER2,
    EXPR_NODEFUNC_MAP,
    EXPR_NODEFUNC_GEOMEDIA,
    EXPR_NODEFUNC_ARMEDIA,
    EXPR_NODEFUNC_POWMEDIA,
    EXPR_NODEFUNC_CVAL,
    EXPR_NODEFUNC_FIRSTQUARTILE,
    EXPR_NODEFUNC_MEDIANA,
    EXPR_NODEFUNC_THIRDQUARTILE,
    //
    EXPR_NODEFUNC_IF,
    EXPR_NODEFUNC_SELECT,
    EXPR_NODEFUNC_EQUAL,
    EXPR_NODEFUNC_ABOVE,
    EXPR_NODEFUNC_BELOW,
    EXPR_NODEFUNC_AVG,
    EXPR_NODEFUNC_CLIP,
    EXPR_NODEFUNC_CLAMP,
    EXPR_NODEFUNC_PNTCHANGE,
    EXPR_NODEFUNC_POLY,
    EXPR_NODEFUNC_AND,
    EXPR_NODEFUNC_OR,
    EXPR_NODEFUNC_NOT,
    EXPR_NODEFUNC_FOR,
    EXPR_NODEFUNC_MANY
    };

/* Forward declarations */
typedef struct _exprFunc exprFunc;
typedef struct _exprVal exprVal;

/* Expression object */
struct _exprObj
    {
    struct _exprFuncList *flist; /* Functions */
    struct _exprValList *vlist; /* Variables */
    struct _exprValList *clist; /* Constants */
    struct _exprNode *headnode; /* Head parsed node */

    exprBreakFuncType breakerfunc; /* Break function type */

    void *userdata; /* User data, can be any 32 bit value */
    int parsedgood; /* non-zero if successfully parsed */
    int parsedbad; /* non-zero if parsed but unsuccessful */
    int breakcount; /* how often to check the breaker function */
    int breakcur; /* do we check the breaker function yet */
    int starterr; /* start position of an error */
    int enderr; /* end position of an error */
    };

/* Object for a function */
struct _exprFunc
    {
    char *fname; /* Name of the function */
    exprFuncType fptr; /* Function pointer */
    int min, max; /* Min and max args for the function. */
    int refmin, refmax; /* Min and max ref. variables for the function */
    int type; /* Function node type.  exprEvalNOde solves the function */

    struct _exprFunc *next; /* For linked list */
    };

/* Function list object */
struct _exprFuncList
    {
    struct _exprFunc *head;
    };

/* Object for values */
struct _exprVal
    {
    char *vname; /* Name of the value */
    EXPRTYPE vval; /* Value of the value */
    EXPRTYPE *vptr; /* Pointer to a value.  Used only if not NULL */

    struct _exprVal *next; /* For linked list */
    };

/* Value list */
struct _exprValList
    {
    struct _exprVal *head;
    };

/* Expression node type */
struct _exprNode
    {
    int type; /* Node type */

    union _data /* Union of info for various types */
        {
        struct _oper
            {
            struct _exprNode *nodes; /* Operation arguments */
            int nodecount; /* Number of arguments */
            } oper;

        struct _variable
            {
            EXPRTYPE *vaddr; /* Used if EXPR_FAST_VAR_ACCESS defined */
            } variable;

        struct _value
            {
            EXPRTYPE value; /* Value if type is value */
            } value;

        struct _assign /* Assignment struct */
            {
            EXPRTYPE *vaddr; /* Used if EXPR_FAST_VAR_ACCESS defined */
            struct _exprNode *node; /* Node to evaluate */
            } assign;

        struct _function
            {
            exprFuncType fptr; /* Function pointer */
            struct _exprNode *nodes; /* Array of argument nodes */
            int nodecount; /* Number of argument nodes */
            EXPRTYPE **refs; /* Reference variables */
            int refcount; /* Number of variable references (not a reference counter) */
            int type; /* Type of function for exprEvalNode if fptr is NULL */
            } function;
        } data;
    };



/* Functions for function lists */
int exprFuncListAddType(exprFuncList *flist, char *name, int type, int min, int max, int refmin, int refmax);
int exprFuncListGet(exprFuncList *flist, char *name, exprFuncType *ptr, int *type, int *min, int *max, int *refmin, int *refmax);


#ifdef __cplusplus
}
#endif

#endif /* __BAVII_EXPRPRIV_H */

