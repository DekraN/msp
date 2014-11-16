/*!////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
//!////////////////////////////////////////////////////////////////////////////////////////////////////////////// ///
/*!________________________________________________________________________________________________________________*/
///                         msp v1.40 --- by Marco Chiarelli aka DekraN aka Wesker013                  	  		  ///
///        							LAST UPDATE: 16:00 of 16/11/2014 by Myself. 								  ///
/// 	This program is protected by Creative Commons CC BY-SA 2.0 License. For more informations contact me. 	  ///
///                     You can contact me at: marco_chiarelli@yahoo.it or on the secundary mail:                 ///
/// marcochiarelli.nextgenlab@gmail.com in order to report a bug or simply for sending me an advice that could be ///
///                        useful or could improve the speed or optimize my MSCenv System.                        ///
//!-------------------------------------------------------------------------------------------------------------- ///
//!______________________________________________________________________________________________________________ ///
//!-------------------------------------------------------------------------------------------------------------- ///
///                contact me at: marco_chiarelli@yahoo.it or marcochiarelli.nextgenlab@gmail.com                 ///
///        I'll be glad to fix your scripts or simply to take away your doubts about the program                  ///
//!-------------------------------------------------------------------------------------------------------------- ///
//!-------------------------------------------------------------------------------------------------------------- ///
/// Thanks to giggikr: http://forum.html.it/forum/showthread/t-1374455.html for his own function, cambiabase,     ///
///                         which I renamed, modified and adapted to this program. Thanks to:                     ///
/// http://elite.polito.it/files/courses/12BHD/progr/Esercizi-C-v2_01.pdf for some of their scripts.              ///
/// Thanks to Bibek Subedi, for his invertMatrix function, which I renamed, modified and adapted to this program. ///
/// Link Source: http://programming-technique.blogspot.it/2011/09/numerical-methods-inverse-of-nxn-matrix.html    ///
/// Thanks to Paul Bourke:https://www.cs.rochester.edu/u/brown/Crypto/assts/projects/adj.html for his CoFactor fnc///
/// 								 which I modified and adapted to this program.								  ///
/// 		Thanks to my University friends: Dino Sbarro, Gabriele Accarino and Giampiero d'Autilia				  ///
/// 			for inspiring me the implementation of some functions and for the everyday support. 			  ///
/// Thanks to W. Cochran  wcochran@vancouver.wsu.edu for his Strassen Algorithm Implementation, which I renamed,  ///
/// adapted and modified to this program. Thanks also to: Computer Science Division | EECS at UC Berkeley for     ///
/// some notions about Matrix Multiplication Optimizations Techniques: www.cs.berkeley.edu/~knight/cs267/hw1.html ///
/// Massive thanks to Brian Allen Vanderburg II for his fabulous C parser and inline functions solver, EXPREVAL,  ///
/// which elegantly gave in theory infinite functionalities and potential to my program. That's the project link  ///
/// with Online Documentation: http://expreval.sourceforge.net/ Thanks to: http://www.cprogramming.com/tips/ and  ///
///             http://stackoverflow.com/questions/599365/what-is-your-favorite-c-programming-trick               ///
///                       http://stackoverflow.com/questions/132241/hidden-features-of-c                          ///
///    that are some websites in which I found a lot of useful C tips and tricks, and they were an important      ///
///  checkpoint for resources retrieving in order to speed-up and optimize my program. Still greatly thanks to    ///
///    Bibek Subedi for his website: http://www.programming-techniques.com/ which put in front of my eyes a new   ///
/// world of C programming. I also recently renewed the program code by improving a lot of his C tricks and tips. ///
/// For example, the upper-triangular Matrixes conversion, which was useful to enhance some functions like det(), ///
///      							sgeqsolver ExprEval inline command, etc.									  ///
///                 Greatly thanks to Daniel Veillard for his fabulous XML Parser, LIBXML2.                       ///
///  Greatly thanks to vict85 of matematicamente.it Network, for informing me about the benefits of using   	  ///
///   generally a single reference for the Matrix Type, like LAPACK and the other Numeric Calculus Environments.  ///
/// 		Thanks to Francesco Palma for reporting me some bugs, and finally, massive thanks to my				  ///
/// Informatic Fundaments Teacher, Mario Alessandro Bochicchio, which gave me a lot of C advices and some general ///
///            tricks and tips, that enlarged my professional informatic horizonts. That's all...                 ///
/*!________________________________________________________________________________________________________________*/
//!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/


#include "dutils.h" // DA RENDERE VISIBILE SIA AL COMPILATORE CHE AL LINKER

struct program * suite;

const struct prog_constants suite_c =
{
    {
    	"Fibonacci",
    	"Factorial",
		"Even SemiFactorial",
		"Odd SemiFactorial",
		"Every"
    },
    {
        "Real Numbers",
        "Complex Numbers",
        "Quaternions",
        "Octonions",
        "Sedenions"
    },
    {
        {
            REAL_UNIT_NAME
        },
        {
            REAL_UNIT_NAME,
            "i"
        },
        {
            REAL_UNIT_NAME,
            "i",
            "j",
            "k"
        },
        {
            REAL_UNIT_NAME,
            "e1",
            "e2",
            "e3",
            "e4",
            "e5",
            "e6",
            "e7"
        },
        {
            REAL_UNIT_NAME,
            "e1",
            "e2",
            "e3",
            "e4",
            "e5",
            "e6",
            "e7",
            "e8",
            "e9",
            "e10",
            "e11",
            "e12",
            "e13",
            "e14",
            "e15"
        }
    }
};

__MSSHELL_WRAPPER_ static void _MS__private __system changeAlgebra(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system switchLogMode(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system progInfo(const sel_typ argc, char ** argv);

sprog main_menu[MAX_PROGRAMMI] =
{
    [MAIN_BCALCULATOR] =
    {
        "Basic Calculator",
        CMD_BCALC,
        USAGE_BCALC,
        basicCalculator,
        BY_USER,
        CHILD
    },
    [MAIN_ADVANCEDCALCULATOR] =
    {
        "Advanced Calculator",
        CMD_ADVCALC,
        USAGE_ADVCALC,
        calcolatoreAvanzato,
        AUTOMATIC,
        FATHER
    },
    [MAIN_ALGEBRAOPERATIONS] =
    {
        "Linear Algebra Operations",
        CMD_LINEARALGEBRAOPERATIONS,
        USAGE_LINEARALGEBRAOPERATIONS,
        algebraOperations,
        AUTOMATIC,
        FATHER
    },
    #ifdef ALLOW_MSSMANAGER
	    [MAIN_MSSMANAGER] =
	    {
	        "PROGRAM Scripts Manager",
	        CMD_MSSMANAGER,
	        USAGE_MSSMANAGER,
	        mssManager,
	        AUTOMATIC,
	        FATHER
	    },
    #endif
    [MAIN_CHANGEALGEBRA] =
    {
        "Change ALGEBRA",
        CMD_CHANGEALGEBRA,
        USAGE_CHANGEALGEBRA,
        changeAlgebra,
        AUTOMATIC,
        CHILD
    },
    [MAIN_SWITCHLOGMODE] =
    {
    	"Switch LOG MODE",
    	CMD_SWITCHLOGMODE,
    	USAGE_SWITCHLOGMODE,
    	switchLogMode,
    	AUTOMATIC,
    	CHILD
    },
    [MAIN_PROGINFO] =
    {
    	"Program INFORMATIONS",
    	CMD_PROGINFO,
    	USAGE_PROGINFO,
    	progInfo
	}
};

__MSSHELL_WRAPPER_ static void _MS__private __system changeAlgebra(const sel_typ argc, char ** argv)
{

    dim_typ i;
    
    if(argc)
    {
        ityp tmp = 0.00;
        if((!parse(argv[0], &tmp)) || tmp != (i = (dim_typ)tmp) || i < MIN_ALGEBRA || i > MAX_ALGEBRA)
        {
            printErr(1, "Invalid inserted Value: not correspondent to any Algebra Identifier");
            printUsage(&main_menu[MAIN_CHANGEALGEBRA]);
            return;
        }
    }
    else if((i = selectListItem(_MAX_ALGEBRA, _MAX_ALGEBRA > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
            "Select Algebra Identifier in which to perform Algebra Operations", suite_c.algebra_elements_names)) == _MAX_ALGEBRA) return;

    access(algebra) = i;
    printf2("%s Algebra has been correctly selected.\n\n", suite_c.algebra_elements_names[i]);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system switchLogMode(const sel_typ argc, char ** argv)
{
	if((access(log_mode) = !access(log_mode)))
	{
		if(!(access(log_pntr) = checkForFHErrors(MSPLOG_PATH, "w")))
			access(log_mode) = false;
	}
	else
		fclose(access(log_pntr));
	#ifdef WINOS
		SetExitButtonState(!access(log_mode));
	#endif
	updInfo();
	printf2("\nLOG MODE has been currently %s.\n\n", access(log_mode) ? "Enabled":"Disabled");
}

__MSSHELL_WRAPPER_ static void _MS__private __system progInfo(const sel_typ argc, char ** argv)
{
    printf2(PROG__NAME" V"PROG__VERSION" by ");
    printf2(PROG__AUTHOR);
    printf2(".\n");
    PRINTL();
    if(argc == WITH_DESCRIPTION)
    {
		printf2("\nThis software is licensed under Creative Commons CC BY-SA 2.0\n");
		printf2("contact me for more informations: marco_chiarelli@yahoo.it\n");
		printf2("or marcochiarelli.nextgenlab@gmail.com. For further informations about\n");
		printf2("ExprEval or MSCenv mathematical set of functions, see basic_calculator.txt\n");
		printf2("or visit the project page at: https://sourceforge.net/projects/msp/\n");
		printf2("NOTE: If you're using this program, then you probably should have downloaded it\n");
		printf2("from the official site or in bundle with mathSuite v6.80 version\n");
		printf2("If not, then you are supposed to contact me ASAP.");
		printf2(EXIT_MESSAGE);
    }
    return;
}

int main(int argc, char **argv)
{
    suite = malloc(sizeof(struct program));
    errMem(suite, HEAPALLOC_ERROR);

    // INITIALIZING suite vars

    access(exprVars) = INIT_VARLIST;

    #pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
    for(dim_typ i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
    {
    	access(sysMem)[i].memoizer = NULL;
    	access(sysMem)[i].current_max_index = 0;
    }

    access(func_list) = INIT_FUNCLIST;
    access(const_list) = INIT_CONSTLIST;

    access(random_seed) = INITIALIZING_RANDOM_SEED;

    __pmode__ = PROGRAM_BUSY;

	access(log_pntr) = NULL;
    access(exitHandle) = INVALID_EXITHANDLE;
    access(log_mode) = INIT_LOGMODE;
    access(sigresult) = INIT_SIGRESULT;
    access(mss) = INIT_MSS;

    access(mode) = PROGRAM_BUSY;
    	
    randomize;
    access(random_seed) = starting_random_seed;

    dim_typ i;

    exprValListFree(access(const_list));
    exprFuncListFree(access(func_list));

    access(func_list) = INIT_FUNCLIST;
    access(const_list) = INIT_CONSTLIST;
    /* Set error buffer */

    int err;
    jmp_buf jumper;

    err = setjmp(jumper);

    if(err)
        {
        /* Free stuff */

        if(access(func_list))
            exprFuncListFree(access(func_list));

        if(err != -1)
            printf2("Error: %d\n", err);

        return EXPREVAL_ERROR;
        }

    err = exprFuncListCreate(&(access(func_list)));
    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Func List Creation Error\n");
        longjmp(jumper, 1);
    }

    /* Init funclist */
    err = exprFuncListInit(access(func_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Error initializing internal functions\n");
        longjmp(jumper, err);
    }

    /* Create constlist */
    err = exprValListCreate(&access(const_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Const List Creation Error\n");
        longjmp(jumper, err);
    }

    /* Init constlist */
    err = exprValListInit(access(const_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Error initializing internal constants\n");
        longjmp(jumper, err);
    }
    
    exprType * const tmp = malloc(sizeof(exprType));
    errMem(tmp, HEAPALLOC_ERROR);

    tmp->var_list = INIT_VARLIST;
    ///    tmp->const_list = INIT_CONSTLIST;
    tmp->e_ANS = INIT_ANS;
    tmp->global_var = DEFAULT_GLOBALVAL;

    /* Create varlist */
    err = exprValListCreate(&(tmp->var_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Var List Creation Error\n");
        longjmp(jumper, err);
    }

    /* Init variable list */
    // err = exprValListAddAddress((*vlist), "global", &(suite.exprVars.global_var));
    err = exprValListAddAddress(tmp->var_list, "global", &(tmp->global_var));
    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Var List Init Error\n");
        longjmp(jumper, err);
    }

    // err = exprValListAdd((*vlist), DEFAULT_ENVS_ANSVALNAME, 0.0);
    err = exprValListAdd(tmp->var_list, DEFAULT_ENVS_ANSVALNAME, 0.00);
    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Error adding variable \'%s\'\n", DEFAULT_ENVS_ANSVALNAME);
        longjmp(jumper, err);
    }

    // exprValListGetAddress((*vlist), DEFAULT_ENVS_ANSVALNAME, &(suite.exprVars.e_res));
    exprValListGetAddress(tmp->var_list, DEFAULT_ENVS_ANSVALNAME, &(tmp->e_ANS));
    if(tmp->e_ANS == NULL)
    {
        printf2("Unable to get address of \'%s\'\n", DEFAULT_ENVS_ANSVALNAME);
        longjmp(jumper, EXPR_ERROR_UNKNOWN);
    }

    access(exprVars) = tmp;
    char * extension = NULL;

    if(argv[1] && strlen(argv[1]) > 1 && (extension= strrchr(argv[1], '.')+1) && !strcmp(extension, DEFAULT_SCRIPTFILES_EXTENSION))
        _execScriptFiles(argv[1]);
    else if(argc > 2)
        _handleCmdLine(argc-2, &argv[2]);

    progInfo(WITHOUT_DESCRIPTION, NULL);

    updInfo();
    operationsGroupMenu(MAX_PROGRAMMI, main_menu, NULL_CHAR, BY_CHARS);

    prepareToExit();
    #ifdef WINOS
    	(void) system("PAUSE"); // if working on Windows Environment...
	#endif

    return NOERROR_EXIT;
}

