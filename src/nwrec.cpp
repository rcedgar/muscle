/***
Needleman-Wunch recursions

Notation: i,j are prefix lengths so are in
ranges i = [0,|A|] and j = [0,|B|].

Profile positions are in ranges [0,|A|-1]
and [0,|B|-1] so prefix length i corresponds
to position (i-1) in the profile, and similarly
for j.

Terminal gap scoring
--------------------
Terminal gaps are scored as with open [close]
penalties only at the left [right] terminal,
as follows:

      0  i
	  |  |
	A XXXXX...
	B ---XX...

          i |A|-1
          |  |
	A ...XXXXX
	B ...XX---

In these examples, open / close penalty at position
i is  included, but close / open penalty at |A|-1 /
0 is not included.

This is implemented by setting the open [close] 
penalty to zero in the first [last] position of
each profile.

Consider adding a column to a sub-alignment. After the
column is added, there are i letters from A and j letters
from B.

The column starts a left-terminal gap if:
	Delete with i=1, j=0 or
	Insert with i=0, j=1.

The column ends a left-terminal gap if:
	Match following Delete with j=1, or
	Match following Insert with i=1.

The column starts a right-terminal gap if:
	Delete following a Match and i=|A|, or
	Insert following a Match and j=|B|.

The column ends a right-terminal gap if:
	Match with i=|A|, j=|B| following Delete or Insert.
	
RECURSION RELATIONS
===================

         i-1
          |
DD	A ..X X
	B ..- -

MD	A ..X X
	B ..X -

D(i,j) = max
			D(i-1,j) + e
			M(i-1,j) + goA(i-1)
Valid for:
	i = [1,|A|-1]
	j = [1,|B|]

I(i,j) By symmetry with D(i,j).

       i-2
        | i-1
		| |
MM	A ..X X
	B ..X X

DM	A ..X X
	B ..- X

IM  A ..- X
	B ..X X
	    | |
		| j-1
	   j-2

M(i,j) = L(i-1,j-1) + max
			M(i-1,j-1)
			D(i-1,j-1) + gcA(i-2)
			I(i-1,j-1) + gcB(j-2)
Valid for:
	i = [2,|A|]
	j = [2,|B|]

Equivalently:

M(i+1,j+1) = L(i,j) + max
			M(i,j)
			D(i,j) + gcA(i-1)
			I(i,j) + gcB(j-1)

Valid for:
	i = [1,|A|-1]
	j = [1,|B|-1]

Boundary conditions
===================

A XXXX
B ----
	D(0,0) = -infinity

	D(i,0) = ie
		i = [1,|A|]

	D(0,j) = -infinity
		j = [0,|B|]

I(0,0), I(0,j) and I(i,0) by symmetry with D.

	M(0,0) = 0
	M(i,0) = -infinity, i > 0
	M(0,j) = -infinity, j > 0

A X
B -
	D(1,0) = e
	D(1,j) = -infinity, j = [1,|B|]
		(assuming no I-D allowed).

	D(0,1) = -infinity
	D(1,1) = -infinity
	D(i,1) = max.
***/
