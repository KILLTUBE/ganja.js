/*
    Table will construct a function like this:

    var arg_0 = 0;
    var arg_1 = 0;
    var arg_2 = 0;
    for (var arg_0 = 0; arg_0 <= 3; arg_0++) {
            for (var arg_1 = 0; arg_1 <= 3; arg_1++) {
                    for (var arg_2 = 0; arg_2 <= 3; arg_2++) {
                            cb(arg_0, arg_1, arg_2);
                    } // arg_2 up to 3
            } // arg_1 up to 3
    } // arg_0 up to 3
*/

/**
 * 
 * @param {number[]} def 
 * @param {()=>void} cb 
 */

export function Table(def, cb) {
    var func = "";
    var i;
    for (i=0; i<def.length; i++) {
        func += "var arg_"+i+" = 0;\n"
    }
    for (i=0; i<def.length; i++) {
        var arg = "arg_" + i;
        var n = def[i];
        var tabs = '\t'.repeat(i);
        func += `${tabs}for (var ${arg} = 0; ${arg} <= ${n}; ${arg}++) {\n`;
    }
    func += '\t'.repeat(def.length);
    func += "cb(";
    for (i in def) {
        var arg = "arg_" + i;
        func += arg;
        if (i < def.length - 1)
            func += ", ";
    }
    func += ");\n";
    for (i=def.length - 1; i>=0; i--) {
        var arg = "arg_" + i;
        var n = def[i];
        var tabs = '\t'.repeat(i);
        func += `${tabs}} // ${arg} up to ${n}\n`;
    }
    //console.log(func);
    var realFunc = new Function("cb", func);
    realFunc(cb)
}

//import {execSync} from 'child_process';

Table([3,3,3], function(i, j, k) {
    // That's cause from 6D Ganja switches to a graded implementation by default
    // (which does no longer store all coefficients of a multivector, only for the used grades)
    // (which btw uses simplify_bits - the bitmask alternative to simplify for its calculations)
    // (and implies that in the graded generator, you no longer get to pick the basis order - its always the grade/lexical ordering)
    // Can be prevented with graded:false option

    if (i + j + k <= 6) {
        var cmd = `node Algebra.js ${i} ${j} ${k} > tests/simplify/${i}_${j}_${k}.js`;
        console.log(cmd)
    }
    //execSync(cmd);
})
