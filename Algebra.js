import Algebra from "./src/Algebra.js";


/*
    process.argv [
        'C:\\node15\\node.exe',
        'D:\\web\\ganja.js\\tests\\simplify_gen.js',
        '3',
        '0',
        '1'
    ]
*/
var p = Number(process.argv[2]);
var q = Number(process.argv[3]);
var r = Number(process.argv[4]);
Algebra(p, q, r);
