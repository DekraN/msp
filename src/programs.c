// programs.c 10/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h" // DA RENDERE VISIBILE       SIA AL COMPILATORE CHE AL LINKER
// #include "ExprEval/exprincl.h" // In order To use Redefined MATH_ constants


__MSSHELL_WRAPPER_ void basicCalculator(const sel_typ argc, char ** argv)
{
	requires(argc ? argv[0] : NULL, argc ? NULL : "Enter an Expression", "Result is", PARSER_SHOWRESULT | PARSER_SHOWVARLIST | PARSER_SHOWDIFFTIME | PARSER_SAVERESULT);
    return ;
}

__MSSHELL_WRAPPER_ void __apnt calcolatoreAvanzato(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_ADVCALC_PROGS, adv_calc, main_menu[MAIN_ADVANCEDCALCULATOR].name, MAX_ADVCALC_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return ;
}

__MSSHELL_WRAPPER_ void __apnt mssManager(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_MSSMANAGER_PROGS, mss_manager, main_menu[MAIN_MSSMANAGER].name, MAX_MSSMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return ;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ void _MS__private __system __export operationsGroupMenu(dim_typ dim, sprog programs[static dim], const char access_point[static INFO_STRING], bool list)
{
    for( ;; )
    {

        printf2("\nSelect desired Program:\n");
        PRINTL();
        // PRINTN();

        dim_typ i;
        dim_typ tmp;

        // tmp = PROGRAM_BUSY;

        CLEARBUFFER();

        if(list == BY_NUMBERS)
        {
            for(i=0; i<dim; ++i)
                printf2("- %hu: %s;\n", i, programs[i].name);

            printf2("- %hu: Clear SCREEN;\n", i);
            printf2("- %hu: PROGRAM Informations;\n", i+1);
            printf2("- %hu: Exit from PROGRAM.\n\n", i+2);
            PRINTL();


            ityp tmp2;

            while(isNullVal((tmp2 = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS))) || tmp2 != (tmp = (dim_typ)tmp2) || tmp < 0 || tmp > dim+2)
                printErr(1, "Invalid PROGRAM Mode");

        }
        else
        {
            // printf2(programs[MAIN_ALGEBRAOPERATIONS].name);
            for(i=0; i<dim; ++i)
                printf2("- %c: %s;\n", i+'A', programs[i].name);
            printf2("- %c: Clear SCREEN;\n", i+'A');
            printf2("- %c: Exit from PROGRAM.\n\n", i+1+'A');
            PRINTL();

            sel_typ tmp2;

            do
                tmp2 = toupper(getch());
            while(tmp2 < 'A' || tmp2 > dim+1+'A');

            // tmp = tmp2-65;
            // tmp = (toupper(tmp2)-65);

            tmp = tmp2-'A';
        }

        __pmode__ = tmp;


        CLEARBUFFER();
        PRINT2N();


        if(tmp == dim)
        {
            pulisciSchermo;
            printf2("\nSCREEN has been correctly cleaned.\n\n");
            PRINTL();
            operationsGroupMenu(dim, programs, access_point, list);
            return;
        }

        if(tmp == dim+1) break;
        sprog prog_chosen = programs[tmp];


        // PRINTING PROGRAM NAME
        char str[INFO_STRING];

        strcpy(str, prog_chosen.name);
        toupper_s(str);
        printf2(str);
        PRINTN();

        bool rep_check;

        do
        {
            // printf2("\n_____________________________________________________\n\n");

            CLEARBUFFER();

            prog_chosen.program_function(0, NULL); // RICHIAMA LA FUNZIONE O METODO DEL PROGRAMMA SELEZIONATO
            
            if((rep_check = (!prog_chosen.isFather) && (!prog_chosen.automatic)))
            {
                PRINTL();
                printf2("Press any key to repeat\nor press %c to go Back to Main Menu.\n", EXIT_CHAR);
            }
            CLEARBUFFER();
        }
        while(rep_check && getch() != EXIT_CHAR);
    }
    
    CLEARBUFFER();
    PRINTL();
    return;
}

// New Families Access Points
//
__MSSHELL_WRAPPER_ void __apnt algebraOperations(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_ALGEBRA_OPERATIONS, alg_operations, main_menu[MAIN_ALGEBRAOPERATIONS].name, MAX_ALGEBRA_OPERATIONS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return;
}
