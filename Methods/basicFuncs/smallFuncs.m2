------------------------------ Small functions that will be used ---------------------------------
-- Cross Product
cross = method(Vector,Vector) := (u,v) -> (
    f := new MutableList from {3:0};
    f#0 = u_(1)*v_(2) - u_(2)*v_(1);
    f#1 = u_(2)*v_(0) - u_(0)*v_(2);
    f#2 = u_(0)*v_(1) - u_(1)*v_(0);
    out := vector(toList f);
    return out
    )

-- Pointwise Multiplication
ptwsProd = method(Vector,Vector) := (u,v) -> (
    return diagonalMatrix(entries(u))*v
    )

-- Dot Product
dot = method(Vector,Vector) := (u,v) -> (
    return sum(entries(diagonalMatrix(entries(u))*v))
    )

-- Scalar Multiplication
scalarProd = method(RingElement,Vector) := (c,v) -> (
    n := #entries(v);
    u := vector({n:c});
    return ptwsProd(u,v)
    )

-- Pointwise divide
ptwsDiv = method(Vector,Vector) := (f,g) ->(
    n = #entries(f);
    h = new MutableList from {n:0};
    for i from 0 to n-1 do(
	    h#i = f_(i)/g_(i);
	  );
    return h;
)

-- Return the diagonal entries of a square matrix
diagEntries = method(Matrix) := (M) -> (
    n = rank source M;
    f := new MutableList from {n:0};
    for i from 0 to n-1 do (
	f#i = M_(i,i);
	);
    return toList f;
    )
---------------------------------------------------------------------------------------------------