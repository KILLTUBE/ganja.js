// import("./src/Algebra.js").then(module=>{var Algebra = module.default; Algebra(3,0,1)});

// node Algebra.js 1 0 2 > tests/simplify_1_0_2.js

function equals(a, b) {
    if (a != b)
        console.error(a, "!=", b);
}