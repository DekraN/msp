// mss_manager.c 10/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifdef ALLOW_MSSMANAGER

__MSSHELL_WRAPPER_ static void _MS__private __system handleCmdLine(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system execScriptFiles(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system showUsage(const sel_typ argc, char ** argv);

sprog mss_manager[MAX_MSSMANAGER_PROGS] =
{
    [MSSMANAGER_COMMANDLINE] =
    {
        "CmdLine Prompt Interpreter",
        CMD_CMDLINE,
        USAGE_CMDLINE,
        handleCmdLine,
        BY_USER, /// AUTOMATIC
        CHILD
    },
    [MSSMANAGER_EXECSCRIPTFILES] =
    {
        "Execute Scriptfile ("DEFAULT_SCRIPTFILES_EXTENSION")",
        CMD_EXECMSSFILES,
        USAGE_EXECMSSFILES,
        execScriptFiles,
        AUTOMATIC,
        CHILD
    },
    [MSSMANAGER_SHOWUSAGES] =
    {
        "CmdLine Informations",
        CMD_SHOWUSAGES,
        USAGE_SHOWUSAGES,
        showUsage, // startingWizard,
        AUTOMATIC,
        CHILD
    }
};

#define MAX_ARGS 5

__MSSHELL_WRAPPER_ static void _MS__private __system handleCmdLine(const sel_typ argc, char ** argv)
{
    printf2("Enter one Programs Macro-List Command.\n\n");

    char str[INFO_STRING];

    (void) gets(str);

    if(!str[0])
        return;

    dim_typ i;
    char *cmdtab[MAX_ARGS];

    for(cmdtab[0]=strtok(str,BLANK_STRING),i=0; cmdtab[i] != NULL; cmdtab[++i] = strtok(NULL, BLANK_STRING));

    _handleCmdLine(i, cmdtab);

    return;
}

__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __system execScriptFiles(const sel_typ argc, char ** argv)
{
    char path[MAX_PATH_LENGTH];
    if(argc)
    {
        if(!file_exists(argv[0]))
        {
            printErr(2, "Inserted Path:\n%s\nrefers to non-existent Scriptfile", argv[0]);
            return;
        }
        strcpy(path, argv[0]);
    }
    else
    {
        bool assert;

        printf2("Enter Path of the "DEFAULT_SCRIPTFILES_EXTENSION" you wish to load.\n");
        printf2("or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
        PRINTL();

        while(scanf("%s", path) != 1 || path[0] == SCANFEXIT_CHAR || (assert = !file_exists(path)))
        {
            CLEARBUFFER();
            if(path[0] == SCANFEXIT_CHAR) return;
            if(assert)
            {
                printErr(2, "Inserted Path:\n%s\nrefers to non-existent Scriptfile", path);
                return;
            }
        // mustcreatefile = true;
        }
    }
    if(_execScriptFiles(path))
        printf2("\nScriptfile:\n%s\nhas been correctly executed.\n\n", path);
    else
        printErr(2, "An error during:\n%s\nScriptfile opening process might have occured", path);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system showUsage(const sel_typ argc, char ** argv)
{
    if(argc)
    {
        const sprog * const prog = searchProgram(argv[0]);
        if(prog)
            _showUsage(prog);
        else
            printErr(1, "SubProgram inserted CmdName hasn't been found into Programs Macro-List");
        return;
    }

    printf2("\n\n*** Program MACRO-LIST ***\n");
    SHOWPAUSEMESSAGE();
    PRINTL();

    dim_typ i;

    for(i=0; i<MAX_PROGRAMMI; ++i)
    {
        _showUsage(&main_menu[i]);
        if(catchPause()) return;
    }

    for(i=0; i<MAX_ADVCALC_PROGS; ++i)
    {
        _showUsage(&adv_calc[i]);
        if(catchPause()) return;
    }

    for(i=0; i<MAX_ALGEBRA_OPERATIONS; ++i)
    {
        _showUsage(&alg_operations[i]);
        if(catchPause()) return;
    }

    return;
}

#endif
