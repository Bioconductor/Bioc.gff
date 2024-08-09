#include "Bioc_gff.h"
#include "readGFF.h"
#include "utils.h"

#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
    /* readGFF.c */
    CALLMETHOD_DEF(gff_colnames, 1),
    CALLMETHOD_DEF(read_gff_pragmas, 1),
    CALLMETHOD_DEF(scan_gff, 5),
    CALLMETHOD_DEF(load_gff, 8),
    /* utils.c */
    CALLMETHOD_DEF(CharacterList_pasteCollapse, 2),
    {NULL, NULL, 0}
};

void R_init_Bioc_gff(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
