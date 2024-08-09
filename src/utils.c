#include "utils.h"

/* Utility for collapsing elements of a CharacterList */
/* Assumes that 'x' is a 'list' */
SEXP CharacterList_pasteCollapse(SEXP x, SEXP sep) {
    SEXP ans;
    if (TYPEOF(x) != VECSXP)
        error("CharacterList_collapse: expected a list");
    PROTECT(ans = allocVector(STRSXP, length(x)));
    for (int i = 0; i < length(x); i++)
        SET_STRING_ELT(ans, i, _STRSXP_collapse(VECTOR_ELT(x, i), sep));
    UNPROTECT(1);
    return ans;
}
