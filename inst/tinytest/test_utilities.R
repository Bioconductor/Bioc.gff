library(IRanges)

pasteCollapse <- Bioc.gff:::pasteCollapse

## Test pasteCollapse
cl <- CharacterList(c("a", "b"), "c", c("d", "e", "f"))
expect_equal(pasteCollapse(cl), c("a,b", "c", "d,e,f"))
expect_equal(pasteCollapse(cl, collapse = "|"), c("a|b", "c", "d|e|f"))

## Test error on wrong input type for x
not_cl <- list(c("a", "b"), "c")
expect_error(
    pasteCollapse(not_cl),
    "'x' must be a CharacterList"
)

## Test error on invalid collapse argument
expect_error(
    pasteCollapse(cl, collapse = ",,"),
    "'collapse' must be a single string, with a single character"
)

## Test with empty CharacterList
cl_empty <- CharacterList()
expect_equal(pasteCollapse(cl_empty), character(0))

## Test with CharacterList containing empty elements
cl_with_empty <- CharacterList(c("a", ""), character(0), c("b", "c"))
expect_equal(pasteCollapse(cl_with_empty), c("a,", "", "b,c"))
