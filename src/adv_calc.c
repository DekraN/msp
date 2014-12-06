// adv_calc.c 10/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"


__MSSHELL_WRAPPER_ static void _MS__private secondGradeEquationSolver(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private complexAdd(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private complexMul(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private getDate(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private polynomEval(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private routhTable(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private juryTable(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private simplexMethod(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private newtonDifferenceTables(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private lagrangeInterpolation(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private funcIntegration(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private straightLineFitting(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private parabolicCurveFitting(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private linearSystemsSolver(const sel_typ argc, char ** argv);


sprog adv_calc[MAX_ADVCALC_PROGS] =
{
    [ADVCALC_SECONDGRADEEQUATIONSOLVER] =
    {
        "Second Grade Equations Solver",
        CMD_SECONDGRADEQSOLVER,
        USAGE_SECONDGRADEQSOLVER,
        secondGradeEquationSolver,
        BY_USER,
        CHILD
    },
    [ADVCALC_COMPLEXNUMBERSSUM] =
    {
        "Complex and HyperComplex Numbers Addition",
        CMD_COMPLEXADD,
        USAGE_COMPLEXADD,
        complexAdd,
        BY_USER,
        CHILD
    },
    [ADVCALC_COMPLEXNUMBERSPROD] =
    {
        "Complex and HyperComplex Numbers Multiplication",
        CMD_COMPLEXMUL,
        USAGE_COMPLEXMUL,
        complexMul,
        BY_USER,
        CHILD
    },
    [ADVCALC_POLYNOMEVALUATOR] =
    {
    	"Polynom Evaluator",
    	CMD_POLYNOMEVALUATOR,
    	USAGE_POLYNOMEVALUATOR,
    	polynomEval,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_POLYNOMDEVALUATOR] =
    {
    	"Polynom Derivative Evaluator",
    	CMD_POLYNOMDEVALUATOR,
    	USAGE_POLYNOMDEVALUATOR,
    	polynomEval,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_ROUTHTABLE] =
    {
    	"Routh Table",
    	CMD_ROUTHTABLE,
    	USAGE_ROUTHTABLE,
    	routhTable,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_JURYTABLE] =
    {
    	"Jury Table",
    	CMD_JURYTABLE,
    	USAGE_JURYTABLE,
    	juryTable,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_SIMPLEXMETHOD] =
    {
        "Non-Dual Simplex Method",
        CMD_SIMPLEXMETHOD,
        USAGE_SIMPLEXMETHOD,
        simplexMethod,
        BY_USER,
        CHILD
    },
    [ADVCALC_NEWTONDIFFTABLES] =
    {
        "Newton Difference Tables",
        CMD_NEWTONDIFFTABLES,
        USAGE_NEWTONDIFFTABLES,
        newtonDifferenceTables,
        BY_USER,
        CHILD
    },
    [ADVCALC_LAGRANGEINTERPOLATION] =
    {
        "Lagrange Unequal Interpolation",
        CMD_LAGRANGEINTERPOLATION,
        USAGE_LAGRANGEINTERPOLATION,
        lagrangeInterpolation,
        BY_USER,
        CHILD
    },
    [ADVCALC_FUNCTIONINTEGRATION] =
    {
        "Function Integration",
        CMD_FID,
        USAGE_FID,
        funcIntegration,
        BY_USER,
        CHILD
    },
    [ADVCALC_STRAIGHTLINEFITTING] =
    {
        "Straight Line Fitting",
        CMD_STRAIGHTLINEFITTING,
        USAGE_STRAIGHTLINEFITTING,
        straightLineFitting,
        BY_USER,
        CHILD
    },
    [ADVCALC_PARABOLICCURVEFITTING] =
    {
        "Parabolic Curve Fitting",
        CMD_PARABOLICCURVEFITTING,
        USAGE_PARABOLICCURVEFITTING,
        parabolicCurveFitting,
        BY_USER,
        CHILD
    },
    [ADVCALC_LINEARSYSTEMSSOLVER] =
    {
        "Linear Systems Solver",
        CMD_LINEARSYSTEMSSOLVER,
        USAGE_LINEARSYSTEMSSOLVER,
        linearSystemsSolver,
        BY_USER,
        CHILD
    }
};

__MSSHELL_WRAPPER_ static void _MS__private secondGradeEquationSolver(const sel_typ argc, char ** argv)
{
    ityp *abc = NULL;

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &abc, dim, &dim[COLUMNS])) || dim[ROWS] != 1 || dim[COLUMNS] != MAX_ABSTRACT_DIMENSIONS)
        {
            matrixFree(&abc);
            printUsage(&adv_calc[ADVCALC_SECONDGRADEEQUATIONSOLVER]);
            return;
        }
    }
    else
    {

        printf2("\nEnter COEFFICIENTS: 'a', 'b' e 'c' of SECOND GRADE Equation:\n\"a*(x^2) + b*x + c\"");
        printf2("\nby inserting related inline [1 x 3] Matrix.\n\n");

        PRINTHOWTOBACKMESSAGE();

        if(!insertNMMatrix(&abc, (dim_typ2){1, MAX_ABSTRACT_DIMENSIONS}))
            return;
    }

    ityp root[MAX_DIMENSIONS];

    if(_secondGradeEquationSolver(abc, root))
    {
        printf2("\n1st ROOT = ");
        printf2(OUTPUT_CONVERSION_FORMAT, root[ROOT_X1]);
        printf2(";\n2nd ROOT = ");
        printf2(OUTPUT_CONVERSION_FORMAT, root[ROOT_X2]);
        printf2(".\n\n");
    }
    
    matrixFree(&abc);

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private complexAdd(const sel_typ argc, char ** argv)
{
    ityp *cpx = NULL;
    const sel_typ algebra_units = (!access(algebra)) ? MAX_COMPLEX_UNITS : exp2(access(algebra));

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &cpx, dim, &dim[COLUMNS])) || dim[ROWS] != MAX_DIMENSIONS || dim[COLUMNS] != algebra_units)
        {
            matrixFree(&cpx);
            printUsage(&adv_calc[ADVCALC_COMPLEXNUMBERSSUM]);
            return;
        }
    }
    else
    {

        printf2("\nEnter 2x%hu MATRIX whose ROWS contains respectively\nREAL PART and IMAGINARY PART%s of both two Operands.\n\n", algebra_units, algebra_units > MAX_COMPLEX_UNITS ? "s":NULL_CHAR);

        if(!insertNMMatrix(&cpx, (dim_typ2){MAX_DIMENSIONS, algebra_units}))
            return;
    }

    ityp complexRes[algebra_units];
    
    static void (* const complexAddFunc[_MAX_ALGEBRA])(ityp *restrict, ityp [static MAX_SEDENIONS_UNITS]) =
	{
		_complexAdd,
		_quaternionsAdd,
		_octonionsAdd,
		_sedenionsAdd
	};
		
	complexAddFunc[((access(algebra) == ALGEBRA_COMPLEXNUMBERS || !access(algebra)) ? ALGEBRA_COMPLEXNUMBERS : access(algebra))-1](cpx, complexRes);
    
	printf2("\nRESULT of Operation: (");
	
    dim_typ i;

    PRINT2N();

    for(i=0; i<algebra_units; ++i)
    {
        printf2(OUTPUT_CONVERSION_FORMAT, *(cpx + i));
        printf2("%s%s ", suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") +\n":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*SECOND_NUMBER) + i));
        printf2("%s%s ", suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") is = to:\n":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(OUTPUT_CONVERSION_FORMAT, complexRes[i]);
        printf2("%s%s", suite_c.algebra_imaginary_units_names[algebra_units][i], i ==algebra_units-1 ? ";\n\n":" + ");
    }

    matrixFree(&cpx);

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private complexMul(const sel_typ argc, char ** argv)
{
    ityp *cpx = NULL;
    const fsel_typ algebra_units = (!access(algebra)) ? MAX_DIMENSIONS : exp2(access(algebra));

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &cpx, dim, &dim[COLUMNS])) || dim[ROWS] != MAX_DIMENSIONS || dim[COLUMNS] != MAX_DIMENSIONS)
        {
            matrixFree(&cpx);
            printUsage(&adv_calc[ADVCALC_COMPLEXNUMBERSPROD]);
            return;
        }
    }
    else
    {

        printf2("\nEnter 2x%hu MATRIX whose ROWS contains respectively\nREAL PART and IMAGINARY PART%s of both two Operands.\n\n", algebra_units, algebra_units > MAX_COMPLEX_UNITS ? "s":NULL_CHAR);

        if(!insertNMMatrix(&cpx, (dim_typ2){MAX_DIMENSIONS, algebra_units}))
            return;
    }

    ityp complexRes[algebra_units];
    
    static void (* const complexMulFunc[_MAX_ALGEBRA])(ityp *restrict, ityp [static MAX_SEDENIONS_UNITS]) =
	{
		_complexMul,
		_quaternionsMul,
		_octonionsMul,
		_sedenionsMul
	};
	
	
	complexMulFunc[((access(algebra) == ALGEBRA_COMPLEXNUMBERS || !access(algebra)) ? ALGEBRA_COMPLEXNUMBERS : access(algebra))-1](cpx, complexRes);
    	
	printf2("\nRESULT of Operation: (");

    dim_typ i;
    
    PRINT2N();

    for(i=0; i<algebra_units; ++i)
    {
        printf2(OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*FIRST_NUMBER) + i));
        printf2("%s%s ", suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") *\n":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*SECOND_NUMBER) + i));
        printf2("%s%s ", suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") is = to:\n":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(OUTPUT_CONVERSION_FORMAT, complexRes[i]);
        printf2("%s%s", suite_c.algebra_imaginary_units_names[algebra_units][i], i ==algebra_units-1 ? ";\n\n":" + ");
    }

    matrixFree(&cpx);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private polynomEval(const sel_typ argc, char ** argv)
{
	ityp *table = NULL;
	dim_typ dim[MAX_DIMENSIONS];
	
	if(argc)
    {
        if((!matrixToken(argv[0], &table, dim, &dim[COLUMNS])))
        {
            matrixFree(&table);
            printUsage(&adv_calc[ADVCALC_POLYNOMEVALUATOR]);
            return;
        }
    }
    else
    {
        printf2("\nEnter the Polynom n-dimensioned Row-Matrix.\n\n");
        if(!insertMatrix(table, dim[ROWS], dim[COLUMNS], false))
            return;
    }
    
    ityp scal = 0.00;

    PRINTHOWTOBACKMESSAGE();
    PRINTN();
    
    if(argc > 1)
    {
	    if(!parse(argv[1], &scal))
	    {
	    	matrixFree(&table);
	    	#ifdef WINOS
	    		SetExitButtonState(ENABLED);
	    	#endif
	        return;
	    }
    }
    else if(isNullVal((scal = requires(NULL, "Enter a double floating-point Scalar Number.\n", "Inserted Scalar", PARSER_SHOWRESULT))))
    {
    	matrixFree(&table);
    	#ifdef WINOS
        	SetExitButtonState(ENABLED);
    	#endif
        return;
    }
    
    dim_typ times = __pmode__ == ADVCALC_POLYNOMDEVALUATOR;
    register ityp result;
    char der_order[MINMIN_STRING] = NULL_CHAR;
    
    if(times)
	{
		ityp tmp = 0.00;
		if(argc > MAX_DIMENSIONS)
		{
            if(!parse(argv[1], &tmp))
            {
            	matrixFree(&table);
				#ifdef WINOS
			    	SetExitButtonState(ENABLED);
			    #endif
                return;
            }
            times = tmp;
	    }
	    else if(isNullVal((tmp = requires(NULL, "Enter a non-zero positive integer representative of the Derivative Order.\n", "Inserted Derivative Order:", PARSER_NOSETTINGS))))
        {
            CLEARBUFFER();
            if(!access(exitHandle))
			{
				#ifdef WINOS
			    	SetExitButtonState(ENABLED);
			    #endif	
			    return;
			}
            printErr(33, "Invalid inserted Value.\nMust be a non-zero integer between 1 and %hu", dim[COLUMNS]);
		}
		
		for(dim_typ i=0; i<times; ++i)
			result = deval(table, dim[COLUMNS], scal);
		printf2("The Derivative of the inserted POLYNOM is the POLYNOM: \n");
		printMatrix(stdout, table, dim);
		sprintf(der_order, "%hu-Derivative ", times);
	}
	else
		result = eval(table, dim[COLUMNS], scal);
    
	printf2("\nInserted POLYNOM %sEvaluated in: ", der_order);
	printf2(OUTPUT_CONVERSION_FORMAT, scal);
	printf2(" RESULT is: ");
	printf2(OUTPUT_CONVERSION_FORMAT, result);
	printf2(".\n\n");
    
    matrixFree(&table);
    
    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS
	return;
}

__MSSHELL_WRAPPER_ static void _MS__private routhTable(const sel_typ argc, char ** argv)
{
	ityp *table = NULL;
	dim_typ dim[MAX_DIMENSIONS];
	
	if(argc)
    {
        if((!matrixToken(argv[0], &table, dim, &dim[COLUMNS])) || dim[COLUMNS] <= 2)
        {
            matrixFree(&table);
            printUsage(&adv_calc[ADVCALC_ROUTHTABLE]);
            return;
        }
    }
    else
    {
        printf2("\nEnter the Polynom n>2 dimensioned Row-Matrix.\n\n");
        if(!insertMatrix(table, dim[ROWS], dim[COLUMNS], false))
            return;
        if(dim[COLUMNS] <= 2)
        {
        	matrixFree(&table);
            printUsage(&adv_calc[ADVCALC_ROUTHTABLE]);
            return;
        }
    }
    
    short permanences;
    fsel_typ nullrow = 0; // could not be null-row-ed the first row
    
    if((permanences = _routhTable(&table, dim[COLUMNS], &nullrow)) == ROUTHTABLE_ALLOC_ERROR)
    	printErr(12, "Routh Table Evaluator Dynamic Memory Allocation Problem");
    else
    {
    	printf2("The ROUTH TABLE of the inserted Polynom Matrix is: \n");
    	printMatrix(stdout, table, (dim_typ2){dim[COLUMNS], ((dim_typ)((dim[COLUMNS]*0.5) + 1))});
    	printf2("PERMANENCES: %hu, VARIATIONS: %hu", permanences, dim[COLUMNS]-1-permanences);
    	if(nullrow)
    		printf2("\nIt has been used the AUXILIARY POLYNOM's Derivative on the %huth NULL ROW.\n\n", nullrow+1);
    }
    
    matrixFree(&table);
    
    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS
	return;
}

__MSSHELL_WRAPPER_ static void _MS__private juryTable(const sel_typ argc, char ** argv)
{
	ityp *table = NULL;
	dim_typ dim[MAX_DIMENSIONS];
	
	if(argc)
    {
        if((!matrixToken(argv[0], &table, dim, &dim[COLUMNS])) || dim[COLUMNS] <= 2)
        {
            matrixFree(&table);
            printUsage(&adv_calc[ADVCALC_JURYTABLE]);
            return;
        }
    }
    else
    {
        printf2("\nEnter the Polynom n>2 dimensioned Row-Matrix.\n\n");
        if(!insertMatrix(table, dim[ROWS], dim[COLUMNS], false))
            return;
        if(dim[COLUMNS] <= 2)
        {
        	matrixFree(&table);
            printUsage(&adv_calc[ADVCALC_JURYTABLE]);
            return;
        }
    }
    
    sel_typ result;
    
    if((result = _juryTable(&table, dim[COLUMNS])) == JURYTABLE_ALLOC_ERROR)
    	printErr(12, "Jury Table Evaluator Dynamic Memory Allocation Problem");
    else
    {
    	printf2("\nThe JURY TABLE of the inserted Polynom Matrix is: \n");
    	printMatrix(stdout, table, (dim_typ2){((dim[COLUMNS]-1)<<1)-3,dim[COLUMNS]});
		printf2("JURY Criterion is: %s.\n", result == JURYTABLE_SATISFIED?"SATISFIED":"NOT SATISFIED");
	}
    
    matrixFree(&table);
    
    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS
	return;
}

__MSSHELL_WRAPPER_ static void _MS__private simplexMethod(const sel_typ argc, char ** argv)
{

    sel_typ mode;
    ityp *tableau = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {
        if((mode = strtod(argv[0], NULL)) == MAX_DIMENSIONS) return;
        if(mode < MIN_PROBLEM || mode > MAX_DIMENSIONS)
        {
            printUsage(&adv_calc[ADVCALC_SIMPLEXMETHOD]);
            return;
        }
    }
    else
    {
        printf2("\nSelect Simplex Method's Problem Type:\n");
        printf2("- A for min Problem,\n- B for max Problem;\n");
        printf2("- %c to go Back...\n\n", EXIT_CHAR);

        do if((mode = toupper(getch())) == EXIT_CHAR) return;
        while(mode < 'A' && mode > EXIT_CHAR);

        mode -= 'A';
    }

    if(argc > 1)
    {
        if(!matrixToken(argv[1], &tableau, dim, &dim[COLUMNS]))
        {
            matrixFree(&tableau);
            printUsage(&adv_calc[ADVCALC_SIMPLEXMETHOD]);
            return;
        }
    }
    else
    {
        printf2("Enter the Constraints' Coefficients Matrix, and use the last row\nto insert the Target Function Coefficients.\n");
        printf2("NOTE: In the last row you have to insert an element\nbefore exiting Insert Process, in order to align Matrix Dimensions.\n\n");
        if(!insertMatrix(tableau, dim[ROWS], dim[COLUMNS], false))
            return;
    }

    dim_typ i;
    dim_typ dimc[MAX_DIMENSIONS];
    const dim_typ dimrows_minus1 = dim[ROWS]-1;

    ityp *constraint_types = NULL;

    if(argc > 2)
    {
        if((!matrixToken(argv[2], &constraint_types, dimc, &dimc[COLUMNS])) || dimc[ROWS] != 1 || dimc[COLUMNS] != dim[ROWS])
        {
            matrixFree(&tableau);
            printUsage(&adv_calc[ADVCALC_SIMPLEXMETHOD]);
            return;
        }
    }
    else
    {
        printf2("\nEnter Constraints Types: 0 for <=, non-zero element for >=.\n");
        if(!insertNMMatrix(&constraint_types, (dim_typ2){1,dimrows_minus1}))
        {
            matrixFree(&tableau);
            return;
        }
    }

    ityp *bfs = NULL;
    const dim_typ bfsdims[MAX_DIMENSIONS] =
    {
        1,
        dim[COLUMNS]-1
    };

    if(!matrixAlloc(&bfs, bfsdims))
    {
        matrixFree(&tableau);
        return;
    }

    sel_typ exit_state;
		
    if((exit_state = _simplexMethod(&tableau, &bfs, dim, constraint_types, mode)) == SIMPLEXMETHOD_INFBFS_ERROR)
        printErr(33, "This Problem has a Solution whose limit is Infinite");
    else if(exit_state == SIMPLEXMETHOD_FARBFS_ERROR)
        printErr(33, "No convergence after %hu iterations!", MAX_SIMPLEXMETHOD_ITERATIONS);
    else if(exit_state == SIMPLEXMETHOD_ALLOC_ERROR)
        printErr(12, "Simplex Method Heap Dynamic Memory Allocation Problem");
    else
    {
        printf2("\nRelaxed Problem BFS with Artificial Variables is: ");
        printMatrix(stdout, bfs, (dim_typ2){1,dim[ROWS]+dim[COLUMNS]-2});
    }
    
    matrixFree(&tableau);
    matrixFree(&constraint_types);
    matrixFree(&bfs);

    return;
}


__MSSHELL_WRAPPER_ static void _MS__private newtonDifferenceTables(const sel_typ argc, char ** argv)
{

    ityp x[MAX_NEWTON_DIFFTABLES_DIM];
    ityp y[MAX_NEWTON_DIFFTABLES_DIM][MAX_NEWTON_DIFFTABLES_DIM];

    dim_typ n;
    ityp tmp;

    if(argc)
    {
        if((!parse(argv[0], &tmp)) || tmp != (n = (dim_typ)tmp) || n < MIN_NEWTON_DIFFTABLES_DIM || n > MAX_NEWTON_DIFFTABLES_DIM)
        {
            printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
            return;
        }
    }
    else
    {
        printf2("\nEnter Difference Table DIMENSION.\n");
        PRINTHOWTOBACKMESSAGE();
        
        while(isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted Difference Table DIMENSION is:", PARSER_NOSETTINGS))) || tmp != (n = (dim_typ)tmp) || n < MIN_NEWTON_DIFFTABLES_DIM || n > MAX_NEWTON_DIFFTABLES_DIM)
        {
            CLEARBUFFER();
            if(!access(exitHandle)) return;
            printErr(5, "Invalid inserted Value.\nMust be a non-negative integer between %hu and %hu", MIN_NEWTON_DIFFTABLES_DIM, MAX_NEWTON_DIFFTABLES_DIM);
        }
    }

    dim_typ i;

    if(argc > 1)
        if(argc == n+1)
        {
            char *token = NULL;
			for(i=0; i<n; ++i)
            {

                if((token = strtok(argv[i+1], ",")))
                {
                    if(!parse(token, &x[i]))
                    {
                        printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                        return;
                    }
                }
                else
                {
                    printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                    return;
                }

                if((token = strtok(NULL, ",")))
                {
                    if(!parse(token, &y[i][0]))
                    {
                        printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                        return;
                    }
                }
                else
                {
                    printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                    return;
                }
            }
        }
        else
        {
            printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
            return;
        }
    else
        for(i=0; i<n; ++i)
        {
            printf2("Enter couple No. %hu as expected format:\n[X]\n[Y]\n", i);
            /// scanf("%f %f",&x[i],&y[i][0]); // PAY STRICT ATTENTION TO THIS HANDLE

            while(isNullVal((x[i] = requires(NULL, NULL_CHAR, "Inserted X is:", PARSER_SHOWRESULT))) ||
                isNullVal((y[i][0] = requires(NULL, NULL_CHAR, "Inserted Y is:", PARSER_SHOWRESULT))))
            {
                CLEARBUFFER();
                if(!access(exitHandle)) return;
                printErr(5, "Invalid inserted Value");
            }
        }
		
    newtonDifferenceTable(n, y, FORWARD_DIFFTAB);
    showNewtonDifferenceTable(n, x, y, FORWARD_DIFFTAB);
    newtonDifferenceTable(n, y, BACKWARD_DIFFTAB);
    showNewtonDifferenceTable(n, x, y, BACKWARD_DIFFTAB);

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private lagrangeInterpolation(const sel_typ argc, char ** argv)
{
    ityp tmp;
    dim_typ dim;

    if(argc)
    {
        if((!parse(argv[0], &tmp)) || tmp != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
        {
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
    }
    else
    {
        printf2("\nEnter Data DIMENSION.\n\n");

        while(((tmp = requires(NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT))) || tmp != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
        {
            CLEARBUFFER();
            if(!access(exitHandle)) return;
            printErr(33, "Invalid inserted Data DIMENSION.\nMust be an integer between 1 and %z", USHRT_MAX);
        }
    }

    ityp *xy = NULL;
    dim_typ i, j;

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy);
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
    }
    else
    {

        printf2("\nEnter related Matrix filled with Data you want Interpolation to Process,\n");
        printf2("by putting on each rows the %hu X and Y Values.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){MAX_DIMENSIONS, dim}))
            return;
    }

    ityp xp;

    if(argc > MAX_DIMENSIONS)
    {
        if(!parse(argv[2], &xp))
        {
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
	}
    else
    {
        printf2("Enter X VALUE to find Y one.\n");
        while(isNullVal((xp = requires(NULL, NULL_CHAR, "X VALUE | Y = F(X) is:", PARSER_SHOWRESULT))))
        {
            CLEARBUFFER();
            if(!access(exitHandle))
                return;
            printErr(5, "Invalid inserted Value");
        }
    }

    ityp yp = 0;
    ityp dr, nr;
		
	for(i = 0; i < dim; ++i)
    {

        dr = nr = 1.00;
        for(j = 0; j<dim; ++j)
            if(i!=j)
            {
                nr *= xp - *(xy + (dim*XROW) + j);
                dr *= *(xy + (dim*XROW) + i) - *(xy + (dim*XROW) + j);
            }

        yp += nr/dr*(*(xy + (dim*YROW) + i));
    }

    matrixFree(&xy);

    printf2("\nRequested Y VALUE is: ");
    printf2(OUTPUT_CONVERSION_FORMAT, yp);
    printf2(".\n\n");

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private funcIntegration(const sel_typ argc, char ** argv)
{
    ityp x0, xn;
    ityp h, s;

    dim_typ funcID;
    dim_typ j;

    funcID = selectListItem(MAX_FIDS, MAX_FIDS > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
                            "Select desired Function you want to Integrate", ext_math.funcnames);

    if(funcID == MAX_FIDS) return;

    ityp tmp;


    sel_typ mode;



    if(argc)
    {
        if((mode = strtod(argv[0], NULL)) == MAX_ABSTRACT_DIMENSIONS) return;
        if(mode < 0 || mode > MAX_ABSTRACT_DIMENSIONS)
        {
            printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
            return;
        }
    }
    else
    {
        printf2("\nSelect Defined Integration Calculus Mode:\n");
        printf2("- A for Simpsons' 1/3 method,\n- B for Simpsons' 3/8 method,\n- C per Trapezoidal Method;\n");
        printf2("- D to go Back...\n\n");

        do if((mode = toupper(getch())) == 'D') return;
        while(mode < 'A' && mode > 'D');

        mode -= 'A';
    }

    dim_typ i, n;

    if(argc > 1)
        if(argc == 4)
        {
            if(!parse(argv[1], &x0))
            {
                printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                return;
            }

            if(!parse(argv[2], &xn))
            {
                printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                return;
            }

            if((!parse(argv[3], &tmp)) || tmp != (n = (dim_typ)tmp))
            {
                printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
                return;
            }
        }
        else
        {
            printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
            return;
        }
    else
    {
        printf2("\nEnter INTEGRATION Extremes and Intervals NUMBER as expected format:\n");
        printf2("[x0]\n[xN]\n[No]\n");

        while(isNullVal((x0 = requires(NULL, NULL_CHAR, "First INTEGRATION Extreme is:", PARSER_SHOWRESULT))) ||
            isNullVal((xn = requires(NULL, NULL_CHAR, "Second INTEGRATION Extreme is:", PARSER_SHOWRESULT))) ||
			isNullVal((tmp = requires(NULL, NULL_CHAR, "Il Numero di INTERVALLI inserito e'", PARSER_SHOWRESULT))) || tmp != (n = (dim_typ)tmp) || n < 1 || n > INT_MAX)
        {
            CLEARBUFFER();
            if(!access(exitHandle)) return;
            printErr(33, "Invalid inserted Intervals NUMBER.\nMust be an integer between 1 and %z", INT_MAX);
        }
    }

    ityp result = 0.00;

    h = (xn - x0) / n;

    ityp (* const y)(register ityp) = ext_math.functions[funcID];

    switch(mode)
    {
        case SIMPSON1DIV8_RULE:
            s = y(x0)+y(xn)+4*y(x0+h);
            for(i = 3; i<=n-1; i+=2)
                s += 4*y(x0+i*h) + 2*y(x0+(i-1)*h);

            result = (h/3)*s;
            break;

        case SIMPSON3DIV8_RULE:
        {
            bool flag;

            s = y(x0)+y(xn);

            for(i = 1; i<=n-1;++i)
            {
                for(j=1;j<=n-1;++j)
                    if((flag = i == 3*j))
                        break;
                s += flag ? 2*y(x0+i*h) : 3*y(x0+i*h);
            }

            result = (3*h/8)*s;
            break;
        }

        case TRAPEZOIDAL_RULE:
            s = y(x0) + y(xn);
            for(i = 0; ++i < n; )
                s += 2*y(x0+i*h);
            result = (h*0.5)*s;
            break;
    }
    
    printf2("%s(x) Function Integral Value calculated between: x0 = ", ext_math.funcnames[funcID]);
    printf2(OUTPUT_CONVERSION_FORMAT, x0);
    printf2(" and xN = ");
    printf2(OUTPUT_CONVERSION_FORMAT, xn);
    printf2(",\nwith %hu Intervals NUMBER is: %6.*f\n\n", n, DEFAULT_PRECISION, result);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private straightLineFitting(const sel_typ argc, char ** argv)
{
    ityp tmp;
    dim_typ dim;

    if(argc)
    {
        if((!parse(argv[0], &tmp)) || tmp != (dim = (dim_typ)tmp))
        {
            printUsage(&adv_calc[ADVCALC_STRAIGHTLINEFITTING]);
            return;
        }
    }
    else
    {
        printf2("\nEnter Data DIMENSION.\n\n");
        while(isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT))) || tmp != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
        {
            CLEARBUFFER();
            if(!access(exitHandle)) return;
            printErr(33, "Invalid inserted Data DIMENSION.\nMust be an integer between 1 and %zu", USHRT_MAX);
        }
    }

    ityp *xy = NULL;

    dim_typ i;

    // we must seek for the BACKTRACKING FEATURE

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy);
            printUsage(&adv_calc[ADVCALC_STRAIGHTLINEFITTING]);
            return;
        }
    }
    else
    {
        printf2("\nEnter the Matrix filled with EQUATIONS DATA,\n");
        printf2("by putting on each ROWS the %hu X and Y VALUES.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){MAX_DIMENSIONS, dim}))
            return;
    }

    ityp sum_x, sum_xy, sum_x2, sum_y;

    sum_x = sum_xy = sum_x2 = sum_y = 0.00;
    
    const register dim_typ cache[MAX_DIMENSIONS] =
    {
    	dim*XROW,
    	dim*YROW
    };
    	
	for(i = 0; i < dim; ++i)
    {

        sum_x += *(xy + cache[XROW] + i);
        sum_y += *(xy + cache[YROW] + i);
        sum_xy += *(xy + cache[XROW] + i) * *(xy + cache[YROW] + i);
        sum_x2 += pow(*(xy + cache[XROW] + i), 2);
    }

    matrixFree(&xy);

    const register ityp b = (dim*sum_xy - sum_x*sum_y)/(dim*sum_x2 - pow(sum_x,2));
    const register ityp a = (sum_y - b*sum_x)/dim;

    printf2("\na = ");
    printf2(OUTPUT_CONVERSION_FORMAT, a);
    printf2("; b = ");
    printf2(OUTPUT_CONVERSION_FORMAT, b);
    printf2(".\nThe Equation is ");
    printf2(OUTPUT_CONVERSION_FORMAT, a);
    printf2(" + ");
    printf2(OUTPUT_CONVERSION_FORMAT, b);
    printf2("X.\n\n");

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private parabolicCurveFitting(const sel_typ argc, char ** argv)
{
    ityp tmp;
    dim_typ dim;

    if(argc)
    {
        if((!parse(argv[0], &tmp)) || tmp != (dim = (dim_typ)tmp))
        {
            printUsage(&adv_calc[ADVCALC_PARABOLICCURVEFITTING]);
            return;
        }
    }
    else
    {
        printf2("\nEnter Data DIMENSION.\n\n");
        while(isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT))) || tmp != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
        {
            CLEARBUFFER();
            if(!access(exitHandle)) return;
            printErr(33, "Invalid inserted Data DIMENSIONM.\nMust be an integer between 1 and %z", USHRT_MAX);
        }
    }

    ityp *xy = NULL;


    // we must seek for the BACKTRACKING FEATURE

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy);
            printUsage(&adv_calc[ADVCALC_PARABOLICCURVEFITTING]);
            return;
        }
    }
    else
    {

        printf2("\nEnter the Matrix filled with EQUATIONS DATA,\n");
        printf2("by putting on each Rows the %hu X and Y VALUES.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){MAX_DIMENSIONS, dim}))
            return;
    }

    ityp *matrix = NULL;

    if(!matrixAlloc(&matrix, (dim_typ2){MAX_ABSTRACT_DIMENSIONS, 4}))
    {
        matrixFree(&xy);
        return;
    }
    
    *(matrix) = dim;

    dim_typ i;
    
    const register dim_typ cache[MAX_DIMENSIONS] =
    {
    	dim*XROW,
    	dim*YROW
    };

    #pragma omp parallel for
	for(i = 0; i < dim; ++i)
    {
        *(matrix + 1) = (*(matrix + 4) += *(xy + cache[XROW] + i));
        *(matrix + 3) += *(xy + cache[YROW] + i);
        *(matrix + 2) = (*(matrix + 8) += pow(*(xy + cache[XROW] + i), 2)); 
        *(matrix + 6) = (*(matrix + 9) += pow(*(xy + cache[XROW] + i), 3));
        *(matrix + 10) += pow(*(xy + cache[XROW] + i), 4); 
        *(matrix + 7) += (*(xy + cache[XROW] + i) * *(xy + cache[YROW] * i));
        *(matrix + 11) += (pow(*(xy + cache[XROW] + i), 2) * *(xy + cache[YROW]*i));
    }

    matrixFree(&xy);

    dim_typ j;
    dim_typ k;

    #pragma omp parallel for
	for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
        #pragma omp parallel for
		for(j = 0; j < MAX_ABSTRACT_DIMENSIONS; ++j)
            if(i!=j)
            {
                const ityp ratio = *(matrix + (j<<MAX_DIMENSIONS) + i)/ *(matrix + (i<<MAX_DIMENSIONS) + i);
                for(k = 0; k < 4; ++k)
                    *(matrix + (j<<MAX_DIMENSIONS) + k) -= ratio * *(matrix + (i<<MAX_DIMENSIONS) + k);
            }

    #pragma omp parallel for
	for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
    {
        const ityp a = *(matrix + (i<<MAX_DIMENSIONS) + i);
		for(j = 0; j < 4; ++j)
            *(matrix + (i<<MAX_DIMENSIONS) + j) /= a;
    }

    printf2("\nPARABOLIC CURVE Fitting is:\n");
    PRINTL();

    for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
    {
        printf2("%c => ", i+97);
        printf2(OUTPUT_CONVERSION_FORMAT, *(matrix + (i<<MAX_DIMENSIONS) + 3));
        printf2(";\n");
    }

    matrixFree(&matrix);

    PRINTL();
    PRINTN();

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private linearSystemsSolver(const sel_typ argc, char ** argv)
{
    dim_typ dim[MAX_DIMENSIONS];
    ityp *matrix = NULL;

    if(argc)
    {

        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[COLUMNS] != dim[ROWS]+1)
        {
            matrixFree(&matrix);
            printUsage(&adv_calc[ADVCALC_LINEARSYSTEMSSOLVER]);
            return;
        }
    }
    else
    {
        printf2("\nEnter Complete Matrix correspondent to the Linear System you want the Program to solve.\n\n");
        if((!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], false)) || dim[COLUMNS] != dim[ROWS]+1)
        {
            if(dim[COLUMNS] != dim[ROWS]+1)
                printErr(33, "You have to insert an [R X R+1] Matrix");
            return;
        }
    }

    dim_typ i, j, k;
    ityp a, b;
    
    for(i = 0; i < dim[ROWS]; ++i)
        for(j = 0; j < dim[ROWS]; ++j)
            if(i != j)
            {
                a = *(matrix + (dim[COLUMNS]*j) + i);
                b = *(matrix + (dim[COLUMNS]*i) + i);
                for(k = 0; k < dim[COLUMNS]; ++k)
                    *(matrix + (dim[COLUMNS]*j) + k) -= (a/b) * *(matrix + (dim[COLUMNS]*i) + k);
            }

	#pragma omp parallel for
    for(i = 0; i < dim[ROWS]; ++i)
    {
        const ityp c = *(matrix + (dim[COLUMNS]*i) + i);
        for(j = 0; j < dim[COLUMNS]; ++j)
            *(matrix + (dim[COLUMNS]*i) + j) /= c;
    }

    printf2("Simultaneous Solutions of the given Linear System are:\n\n");
    PRINTL();

    for(i = 0; i < dim[ROWS] ; ++i)
    {
        printf2("%c => ", i+97);
        printf2(OUTPUT_CONVERSION_FORMAT, *(matrix + (dim[COLUMNS]*i) + dim[ROWS]));
        printf2(";\n");
    }

    matrixFree(&matrix);

    PRINTL();
    PRINTN();
    return;
}
