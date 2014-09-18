#include "dutils.h"

__MSSHELL_WRAPPER_ static void _MS__private __system changeAlgebra(const sel_typ argc, char ** argv);

__MSSHELL_WRAPPER_ static void _MS__private __system changeAlgebra(const sel_typ argc, char ** argv)
{

    dim_typ i;
    
    if(argc)
    {
        ityp tmp = 0.00;
        if((!parse(argv[0], &tmp)) || tmp != (i = (dim_typ)tmp) || i < MIN_ALGEBRA || i > MAX_ALGEBRA)
        {
            printErr(1, "Invalid inserted Value: not correspondent to any Algebra Identifier");
            printUsage(&change_settings[SETTINGS_CHANGEALGEBRA]);
            return;
        }
    }
    else if((i = selectListItem(_MAX_ALGEBRA, _MAX_ALGEBRA > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
            "Select Algebra Identifier in which to perform Algebra Operations", suite_c.algebra_elements_names)) == _MAX_ALGEBRA) return;

    access(algebra) = i;
    printf("%s Algebra has been correctly selected.\n\n", suite_c.algebra_elements_names[i]);
    return;
}
