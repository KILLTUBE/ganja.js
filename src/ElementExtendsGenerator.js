import Element_graph from "./Element_graph.js";
import Element_graphGL from "./Element_graphGL.js";
import Element_graphGL2 from "./Element_graphGL2.js";

/**
 * 
 * @param {*} generator 
 * @param {Options} options 
 * @param {*} tot 
 * @param {*} drm 
 * @param {*} counts 
 * @param {*} simplify 
 * @param {*} basis 
 * @param {*} p 
 * @param {*} q 
 * @param {*} r 
 * @param {*} metric 
 * @param {*} mulTable 
 * @param {*} grades 
 * @param {*} gp 
 * @param {*} basisg 
 */

export default function ElementExtendsGenerator(
    generator,
    options,
    tot,
    drm,
    counts,
    simplify,
    basis,
    p,
    q,
    r,
    metric,
    mulTable,
    grades,
    gp,
    basisg
) {
    // Generate a new class for our algebra. It extends the javascript typed arrays (default float32 but can be specified in options).
    var res = class Element extends generator {

        // constructor - create a floating point array with the correct number of coefficients.
        /**
         * 
         * @param {*} [a]
         */
        constructor(a) { super(a); if (this.upgrade) this.upgrade(); return this; }

        // Grade selection. (implemented by parent class).
        /**
         * 
         * @param {*} grade 
         * @param {*} res 
         */
        Grade(grade, res) { res = res || new Element(); return super.Grade(grade, res); }

        // Right and Left divide - Defined on the elements, shortcuts to multiplying with the inverse.
        Div(b, res) { return this.Mul(b.Inverse, res); }
        LDiv(b, res) { return b.Inverse.Mul(this, res); }

        // Taylor exp - for PGA bivectors in 2D and 3D closed form solution is used.
        Exp() {
            if (options.dual) { var f = Math.exp(this.s); return this.map((x, i) => i ? x * f : f); }
            if (r == 1 && tot <= 4 && Math.abs(this[0]) < 1E-9 && !options.over) {
                var u = Math.sqrt(Math.abs(this.Dot(this).s)); if (Math.abs(u) < 1E-5) return this.Add(Element.Scalar(1));
                var v = this.Wedge(this).Scale(-1 / (2 * u));
                var res2 = Element.Add(Element.Sub(Math.cos(u), v.Scale(Math.sin(u))), Element.Div(Element.Mul((Element.Add(Math.sin(u), v.Scale(Math.cos(u)))), this), (Element.Add(u, v))));
                return res2;
            }
            var res = Element.Scalar(1), y = 1, M = this.Scale(1), N = this.Scale(1); for (var x = 1; x < 15; x++) { res = res.Add(M.Scale(1 / y)); M = M.Mul(N); y = y * (x + 1); }; return res;
        }

        // Log - only for up to 3D PGA for now
        Log() {
            if (r != 1 || tot > 4 || options.over) return;
            var b = this.Grade(2), bdb = Element.Dot(b, b);
            if (Math.abs(bdb.s) <= 1E-5) return this.s < 0 ? b.Scale(-1) : b;
            var s = Math.sqrt(-bdb), bwb = Element.Wedge(b, b);
            if (Math.abs(bwb.e0123) <= 1E-5) return b.Scale(Math.atan2(s, this.s) / s);
            var p = bwb.Scale(-1 / (2 * s));
            return Element.Div(Element.Mul(Element.Mul((Element.Add(Math.atan2(s, this.s), Element.Div(p, this.s))), b), (Element.Sub(s, p))), (Element.Mul(s, s)));
        }

        // Helper for efficient inverses. (custom involutions - negates grades in arguments).
        Map() { var res = new Element(); return super.Map(res, ...arguments); }

        // Factories - Make it easy to generate vectors, bivectors, etc when using the functional API. None of the examples use this but
        // users that have used other GA libraries will expect these calls. The Coeff() is used internally when translating algebraic literals.
        static Element() { return new Element([...arguments]); };
        static Coeff() { return (new Element()).Coeff(...arguments); }
        static Scalar(x) { return (new Element()).Coeff(0, x); }
        static Vector() { return (new Element()).nVector(1, ...arguments); }
        static Bivector() { return (new Element()).nVector(2, ...arguments); }
        static Trivector() { return (new Element()).nVector(3, ...arguments); }
        static nVector(n) { return (new Element()).nVector(...arguments); }

        // Static operators. The parser will always translate operators to these static calls so that scalars, vectors, matrices and other non-multivectors can also be handled.
        // The static operators typically handle functions and matrices, calling through to element methods for multivectors. They are intended to be flexible and allow as many
        // types of arguments as possible. If performance is a consideration, one should use the generated element methods instead. (which only accept multivector arguments)
        static toEl(x) { if (x instanceof Function) x = x(); if (!(x instanceof Element)) x = Element.Scalar(x); return x; }

        // Addition and subtraction. Subtraction with only one parameter is negation.
        static Add(a, b, res) {
            // Resolve expressions passed in.
            while (a.call) a = a(); while (b.call) b = b(); if (a.Add && b.Add) return a.Add(b, res);
            // If either is a string, the result is a string.
            if ((typeof a == 'string') || (typeof b == 'string')) return a.toString() + b.toString();
            // If only one is an array, add the other element to each of the elements.
            if ((a instanceof Array && !a.Add) ^ (b instanceof Array && !b.Add)) return (a instanceof Array) ? a.map(x => Element.Add(x, b)) : b.map(x => Element.Add(a, x));
            // If both are equal length arrays, add elements one-by-one
            if ((a instanceof Array) && (b instanceof Array) && a.length == b.length) return a.map((x, xi) => Element.Add(x, b[xi]));
            // If they're both not elements let javascript resolve it.
            if (!(a instanceof Element || b instanceof Element)) return a + b;
            // Here we're left with scalars and multivectors, call through to generated code.
            a = Element.toEl(a); b = Element.toEl(b); return a.Add(b, res);
        }

        static Sub(a, b, res) {
            // Resolve expressions passed in.
            while (a.call) a = a(); while (b && b.call) b = b(); if (a.Sub && b && b.Sub) return a.Sub(b, res);
            // If only one is an array, add the other element to each of the elements.
            if (b && ((a instanceof Array) ^ (b instanceof Array))) return (a instanceof Array) ? a.map(x => Element.Sub(x, b)) : b.map(x => Element.Sub(a, x));
            // If both are equal length arrays, add elements one-by-one
            if (b && (a instanceof Array) && (b instanceof Array) && a.length == b.length) return a.map((x, xi) => Element.Sub(x, b[xi]));
            // Negation
            if (arguments.length == 1) return Element.Mul(a, -1);
            // If none are elements here, let js do it.
            if (!(a instanceof Element || b instanceof Element)) return a - b;
            // Here we're left with scalars and multivectors, call through to generated code.
            a = Element.toEl(a); b = Element.toEl(b); return a.Sub(b, res);
        }

        // The geometric product. (or matrix*matrix, matrix*vector, vector*vector product if called with 1D and 2D arrays)
        static Mul(a, b, res) {
            // Resolve expressions
            while (a.call && !a.length) {
                a = a();
            }
            while (b.call && !b.length) {
                b = b();
            }
            if (a.Mul && b.Mul) {
                return a.Mul(b, res);
            }
            // still functions -> experimental curry style (dont use this.)
            if (a.call && b.call) return (ai, bi) => Element.Mul(a(ai), b(bi));
            // scalar mul.
            if (Number.isFinite(a) && b.Scale) return b.Scale(a); else if (Number.isFinite(b) && a.Scale) return a.Scale(b);
            // Handle matrices and vectors.
            if ((a instanceof Array) && (b instanceof Array)) {
                // vector times vector performs a dot product. (which internally uses the GP on each component)
                if ((!(a[0] instanceof Array) || (a[0] instanceof Element)) && (!(b[0] instanceof Array) || (b[0] instanceof Element))) { var r = tot ? Element.Scalar(0) : 0; a.forEach((x, i) => r = Element.Add(r, Element.Mul(x, b[i]), r)); return r; }
                // Array times vector
                if (!(b[0] instanceof Array)) return a.map((x, i) => Element.Mul(a[i], b));
                // Array times Array
                var r = a.map((x, i) => b[0].map((y, j) => { var r = tot ? Element.Scalar(0) : 0; x.forEach((xa, k) => r = Element.Add(r, Element.Mul(xa, b[k][j]))); return r; }));
                // Return resulting array or scalar if 1 by 1.
                if (r.length == 1 && r[0].length == 1) return r[0][0]; else return r;
            }
            // Only one is an array multiply each of its elements with the other.
            if ((a instanceof Array) ^ (b instanceof Array)) return (a instanceof Array) ? a.map(x => Element.Mul(x, b)) : b.map(x => Element.Mul(a, x));
            // Try js multiplication, else call through to geometric product.
            var r = a * b; if (!isNaN(r)) return r;
            a = Element.toEl(a); b = Element.toEl(b); return a.Mul(b, res);
        }

        // The inner product. (default is left contraction).
        static LDot(a, b, res) {
            // Expressions
            while (a.call) a = a(); while (b.call) b = b(); //if (a.LDot) return a.LDot(b,res);
            // Map elements in array
            if (b instanceof Array && !(a instanceof Array)) return b.map(x => Element.LDot(a, x));
            if (a instanceof Array && !(b instanceof Array)) return a.map(x => Element.LDot(x, b));
            // js if numbers, else contraction product.
            if (!(a instanceof Element || b instanceof Element)) return a * b;
            a = Element.toEl(a); b = Element.toEl(b); return a.LDot(b, res);
        }

        // The symmetric inner product. (default is left contraction).
        static Dot(a, b, res) {
            // Expressions
            while (a.call) a = a(); while (b.call) b = b(); //if (a.LDot) return a.LDot(b,res);
            // js if numbers, else contraction product.
            if (!(a instanceof Element || b instanceof Element)) return a | b;
            a = Element.toEl(a); b = Element.toEl(b); return a.Dot(b, res);
        }

        // The outer product. (Grassman product - no use of metric)
        static Wedge(a, b, res) {
            // Expressions
            while (a.call) a = a(); while (b.call) b = b(); if (a.Wedge) return a.Wedge(Element.toEl(b), res);
            // The outer product of two vectors is a matrix .. internally Mul not Wedge !
            if (a instanceof Array && b instanceof Array) return a.map(xa => b.map(xb => Element.Mul(xa, xb)));
            // js, else generated wedge product.
            if (!(a instanceof Element || b instanceof Element)) return a * b;
            a = Element.toEl(a); b = Element.toEl(b); return a.Wedge(b, res);
        }

        // The regressive product. (Dual of the outer product of the duals).
        static Vee(a, b, res) {
            // Expressions
            while (a.call) a = a(); while (b.call) b = b(); if (a.Vee) return a.Vee(Element.toEl(b), res);
            // js, else generated vee product. (shortcut for dual of wedge of duals)
            if (!(a instanceof Element || b instanceof Element)) return 0;
            a = Element.toEl(a); b = Element.toEl(b); return a.Vee(b, res);
        }

        // The sandwich product. Provided for convenience (>>> operator)
        static sw(a, b) {
            // Expressions
            while (a.call) a = a(); while (b.call) b = b(); if (a.sw) return a.sw(b);
            // Map elements in array
            if (b instanceof Array && !b.Add) return b.map(x => Element.sw(a, x));
            // Call through. no specific generated code for it so just perform the muls.
            a = Element.toEl(a); b = Element.toEl(b); return a.Mul(b).Mul(a.Reverse);
        }

        // Division - scalars or cal through to element method.
        static Div(a, b, res) {
            // Expressions
            while (a.call) a = a(); while (b.call) b = b();
            // For DDG experiments, I'll include a quick cholesky on matrices here. (vector/matrix)
            if ((a instanceof Array) && (b instanceof Array) && (b[0] instanceof Array)) {
                // factor
                var R = b.flat(), i, j, k, sum, i_n, j_n, n = b[0].length, s = new Array(n), x = new Array(n), yi;
                for (i = 0; i < n; i++) {
                    i_n = i * n;
                    for (j = 0; j < i; j++) {
                        j_n = j * n;
                        s[j] = R[i_n + j];
                        for (k = 0; k < j; k++) s[j] -= s[k] * R[j_n + k];
                        if (R[j_n + j] == 0) return null;
                        R[i_n + j] = s[j] / R[j_n + j];
                    }
                    sum = R[i_n + i];
                    for (k = 0; k < i; k++)  sum -= s[k] * R[i_n + k];
                    R[i_n + i] = sum;
                }
                // subst
                for (i = 0; i < n; i++) for (x[i] = a[i], j = 0; j <= i - 1; j++) x[i] -= R[i * n + j] * x[j];
                for (i = n - 1; i >= 0; i--) for (x[i] /= R[i * n + i], j = i + 1; j < n; j++) x[i] -= R[j * n + i] * x[j];
                return x;
            }
            // js or call through to element divide.
            if (!(a instanceof Element || b instanceof Element)) return a / b;
            a = Element.toEl(a);
            if (Number.isFinite(b)) { return a.Scale(1 / b, res); }
            b = Element.toEl(b); return a.Div(b, res);
        }

        // Pow - needs obvious extensions for natural powers. (exponentiation by squaring)
        static Pow(a, b, res) {
            // Expressions
            while (a.call) a = a(); while (b.call) b = b(); if (a.Pow) return a.Pow(b, res);
            // Exponentiation.
            if (a === Math.E && b.Exp) return b.Exp();
            // Squaring
            if (b === 2) return this.Mul(a, a, res);
            // No elements, call through to js
            if (!(a instanceof Element || b instanceof Element)) return a ** b;
            // Inverse
            if (b === -1) return a.Inverse;
            // Call through to element pow.
            a = Element.toEl(a); return a.Pow(b);
        }

        // Handles scalars and calls through to element method.
        static exp(a) {
            // Expressions.
            while (a.call) a = a();
            // If it has an exp callthrough, use it, else call through to math.
            if (a.Exp) return a.Exp();
            return Math.exp(a);
        }

        // Dual, Involute, Reverse, Conjugate, Normalize and length, all direct call through. Conjugate handles matrices.
        static Dual(a) { return Element.toEl(a).Dual; };
        static Involute(a) { return Element.toEl(a).Involute; };
        static Reverse(a) { return Element.toEl(a).Reverse; };
        static Conjugate(a) { if (a.Conjugate) return a.Conjugate; if (a instanceof Array) return a[0].map((c, ci) => a.map((r, ri) => Element.Conjugate(a[ri][ci]))); return Element.toEl(a).Conjugate; }
        static Normalize(a) { return Element.toEl(a).Normalized; };
        static Length(a) { return Element.toEl(a).Length };

        // Comparison operators always use length. Handle expressions, then js or length comparison
        static eq(a, b) { if (!(a instanceof Element) || !(b instanceof Element)) return a == b; while (a.call) a = a(); while (b.call) b = b(); for (var i = 0; i < a.length; i++) if (a[i] != b[i]) return false; return true; }
        static neq(a, b) { if (!(a instanceof Element) || !(b instanceof Element)) return a != b; while (a.call) a = a(); while (b.call) b = b(); for (var i = 0; i < a.length; i++) if (a[i] != b[i]) return true; return false; }
        static lt(a, b) { while (a.call) a = a(); while (b.call) b = b(); return (a instanceof Element ? a.Length : a) < (b instanceof Element ? b.Length : b); }
        static gt(a, b) { while (a.call) a = a(); while (b.call) b = b(); return (a instanceof Element ? a.Length : a) > (b instanceof Element ? b.Length : b); }
        static lte(a, b) { while (a.call) a = a(); while (b.call) b = b(); return (a instanceof Element ? a.Length : a) <= (b instanceof Element ? b.Length : b); }
        static gte(a, b) { while (a.call) a = a(); while (b.call) b = b(); return (a instanceof Element ? a.Length : a) >= (b instanceof Element ? b.Length : b); }

        // Debug output and printing multivectors.
        static describe(x) { if (x === true) console.log(`Basis\n${basis}\nMetric\n${metric.slice(1, 1 + tot)}\nCayley\n${mulTable.map(x => (x.map(x => ('           ' + x).slice(-2 - tot)))).join('\n')}\nMatrix Form:\n` + gp.map(x => x.map(x => x.match(/(-*b\[\d+\])/)).map(x => x && ((x[1].match(/-/) || ' ') + String.fromCharCode(65 + 1 * x[1].match(/\d+/))) || ' 0')).join('\n')); return { basis: basisg || basis, metric, mulTable, matrix: gp.map(x => x.map(x => x.replace(/\*this\[.+\]/, '').replace(/b\[(\d+)\]/, (a, x) => (metric[x] == -1 || metric[x] == 0 && grades[x] > 1 && (-1) ** grades[x] == (metric[basis.indexOf(basis[x].replace('0', ''))] || (-1) ** grades[x]) ? '-' : '') + basis[x]).replace('--', ''))) } }

        // Direct sum of algebras - experimental
        static sum(B) {
            var A = Element;
            // Get the multiplication tabe and basis.
            var T1 = A.describe().mulTable, T2 = B.describe().mulTable;
            var B1 = A.describe().basis, B2 = B.describe().basis;
            // Get the maximum index of T1, minimum of T2 and rename T2 if needed.
            var max_T1 = B1.filter(x => x.match(/e/)).map(x => x.match(/\d/g)).flat().map(x => x | 0).sort((a, b) => b - a)[0];
            var max_T2 = B2.filter(x => x.match(/e/)).map(x => x.match(/\d/g)).flat().map(x => x | 0).sort((a, b) => b - a)[0];
            var min_T2 = B2.filter(x => x.match(/e/)).map(x => x.match(/\d/g)).flat().map(x => x | 0).sort((a, b) => a - b)[0];
            // remapping ..
            T2 = T2.map(x => x.map(y => y.match(/e/) ? y.replace(/(\d)/g, (x) => (x | 0) + max_T1) : y.replace("1", "e" + (1 + max_T2 + max_T1))));
            B2 = B2.map((y, i) => i == 0 ? y.replace("1", "e" + (1 + max_T2 + max_T1)) : y.replace(/(\d)/g, (x) => (x | 0) + max_T1));
            // Build the new basis and multable..
            var basis = [...B1, ...B2];
            var Cayley = T1.map((x, i) => [...x, ...T2[0].map(x => "0")]).concat(T2.map((x, i) => [...T1[0].map(x => "0"), ...x]))
            // Build the new algebra.
            var grades = [...B1.map(x => x == "1" ? 0 : x.length - 1), ...B2.map((x, i) => i ? x.length - 1 : 0)];
            var a = Algebra({ basis, Cayley, grades, tot: Math.log2(B1.length) + Math.log2(B2.length) })
            // And patch up ..
            a.Scalar = function (x) {
                var res = new a();
                for (var i = 0; i < res.length; i++) res[i] = basis[i] == Cayley[i][i] ? x : 0;
                return res;
            }
            return a;
        }

        // The inline function is a js to js translator that adds operator overloading and algebraic literals.
        // It can be called with a function, a string, or used as a template function.
        static inline(intxt) {
            // If we are called as a template function.
            if (arguments.length > 1 || intxt instanceof Array) {
                var args = [].slice.call(arguments, 1);
                return res.inline(new Function(args.map((x, i) => '_template_' + i).join(), 'return (' + intxt.map((x, i) => (x || '') + (args[i] && ('_template_' + i) || '')).join('') + ')')).apply(res, args);
            }
            // Get the source input text.
            var txt = (intxt instanceof Function) ? intxt.toString() : `function(){return (${intxt})}`;
            // Our tokenizer reads the text token by token and stores it in the tok array (as type/token tuples).
            var tok = [], resi = [], t, tokens = [/^[\s\uFFFF]|^[\u000A\u000D\u2028\u2029]|^\/\/[^\n]*\n|^\/\*[\s\S]*?\*\//g,                 // 0: whitespace/comments
                /^\"\"|^\'\'|^\".*?[^\\]\"|^\'.*?[^\\]\'|^\`[\s\S]*?[^\\]\`/g,                                                                // 1: literal strings
                /^\d+[.]{0,1}\d*[ei][\+\-_]{0,1}\d*|^\.\d+[ei][\+\-_]{0,1}\d*|^e_\d*/g,                                                       // 2: literal numbers in scientific notation (with small hack for i and e_ asciimath)
                /^\d+[.]{0,1}\d*[E][+-]{0,1}\d*|^\.\d+[E][+-]{0,1}\d*|^0x\d+|^\d+[.]{0,1}\d*|^\.\d+|^\(\/.*[^\\]\/\)/g,                       // 3: literal hex, nonsci numbers and regex (surround regex with extra brackets!)
                /^(\.Normalized|\.Length|\.\.\.|>>>=|===|!==|>>>|<<=|>>=|=>|\|\||[<>\+\-\*%&|^\/!\=]=|\*\*|\+\+|\-\-|<<|>>|\&\&|\^\^|^[{}()\[\];.,<>\+\-\*%|&^!~?:=\/]{1})/g,   // 4: punctuator
                /^[A-Za-z0-9_]*/g]                                                                                                            // 5: identifier
            while (txt.length) for (t in tokens) if (resi = txt.match(tokens[t])) { tok.push([t | 0, resi[0]]); txt = txt.slice(resi[0].length); break; } // tokenise
            // Translate algebraic literals. (scientific e-notation to "this.Coeff"
            tok = tok.map(t => (t[0] == 2) ? [2, 'Element.Coeff(' + basis.indexOf((!options.Cayley ? simplify : (x) => x)('e' + t[1].split(/e_|e|i/)[1] || 1).replace('-', '')) + ',' + (simplify(t[1].split(/e_|e|i/)[1] || 1).match('-') ? "-1*" : "") + parseFloat(t[1][0] == 'e' ? 1 : t[1].split(/e_|e|i/)[0]) + ')'] : t);
            // String templates (limited support - needs fundamental changes.).
            tok = tok.map(t => (t[0] == 1 && t[1][0] == '`') ? [1, t[1].replace(/\$\{(.*?)\}/g, a => "${" + Element.inline(a.slice(2, -1)).toString().match(/return \((.*)\)/)[1] + "}")] : t);
            // We support two syntaxes, standard js or if you pass in a text, asciimath.
            var syntax = (intxt instanceof Function) ? [[['.Normalized', 'Normalize', 2], ['.Length', 'Length', 2]], [['~', 'Conjugate', 1], ['!', 'Dual', 1]], [['**', 'Pow', 0, 1]], [['^', 'Wedge'], ['&', 'Vee'], ['<<', 'LDot']], [['*', 'Mul'], ['/', 'Div']], [['|', 'Dot']], [['>>>', 'sw', 0, 1]], [['-', 'Sub'], ['+', 'Add']], [['%', '%']], [['==', 'eq'], ['!=', 'neq'], ['<', 'lt'], ['>', 'gt'], ['<=', 'lte'], ['>=', 'gte']]]
                : [[['pi', 'Math.PI'], ['sin', 'Math.sin']], [['ddot', 'this.Reverse'], ['tilde', 'this.Involute'], ['hat', 'this.Conjugate'], ['bar', 'this.Dual']], [['^', 'Pow', 0, 1]], [['^^', 'Wedge'], ['*', 'LDot']], [['**', 'Mul'], ['/', 'Div']], [['-', 'Sub'], ['+', 'Add']], [['<', 'lt'], ['>', 'gt'], ['<=', 'lte'], ['>=', 'gte']]];
            // For asciimath, some fixed translations apply (like pi->Math.PI) etc ..
            tok = tok.map(t => (t[0] != 5) ? t : [].concat.apply([], syntax).filter(x => x[0] == t[1]).length ? [5, [].concat.apply([], syntax).filter(x => x[0] == t[1])[0][1]] : t);
            // Now the token-stream is translated recursively.
            function translate(tokens) {
                // helpers : first token to the left of x that is not of a type in the skip list.
                var left = (x = ti - 1, skip = [0]) => { while (x >= 0 && ~skip.indexOf(tokens[x][0])) x--; return x; },
                    // first token to the right of x that is not of a type in the skip list.
                    right = (x = ti + 1, skip = [0]) => { while (x < tokens.length && ~skip.indexOf(tokens[x][0])) x++; return x; },
                    // glue from x to y as new type, optionally replace the substring with sub.
                    glue = (x, y, tp = 5, sub) => { tokens.splice(x, y - x + 1, [tp, ...(sub || tokens.slice(x, y + 1))]) },
                    // match O-C pairs. returns the 'matching bracket' position
                    match = (O = "(", C = ")") => { var o = 1, x = ti + 1; while (o) { if (tokens[x][1] == O) o++; if (tokens[x][1] == C) o--; x++; }; return x - 1; };
                // grouping (resolving brackets).
                for (var ti = 0, t, si; t = tokens[ti]; ti++) if (t[1] == "(") glue(ti, si = match(), 6, [[4, "("], ...translate(tokens.slice(ti + 1, si)), [4, ")"]]);
                // [] dot call and new
                for (var ti = 0, t, si; t = tokens[ti]; ti++) {
                    if (t[1] == "[") { glue(ti, si = match("[", "]"), 6, [[4, "["], ...translate(tokens.slice(ti + 1, si)), [4, "]"]]); if (ti) ti--; }    // matching []
                    else if (t[1] == ".") { glue(left(), right()); ti--; }                                                                   // dot operator
                    else if (t[0] == 6 && ti && left() >= 0 && tokens[left()][0] >= 5 && tokens[left()][1] != "return") { glue(left(), ti--) }     // collate ( and [
                    else if (t[1] == 'new') { glue(ti, right()) };                                                                           // collate new keyword
                }
                // ++ and --
                for (var ti = 0, t; t = tokens[ti]; ti++) if (t[1] == "++" || t[1] == "--") glue(left(), ti);
                // unary - and + are handled seperately from syntax ..
                for (var ti = 0, t, si; t = tokens[ti]; ti++)
                    if (t[1] == "-" && (left() < 0 || (tokens[left()] || [4])[0] == 4)) glue(ti, right(), 5, ["Element.Sub(", tokens[right()], ")"]);   // unary minus works on all types.
                    else if (t[1] == "+" && (tokens[left()] || [0])[0] == 4 && (tokens[left()] || [0])[1][0] != ".") glue(ti, ti + 1);                   // unary plus is glued, only on scalars.
                // now process all operators in the syntax list ..
                for (var si = 0, s; s = syntax[si]; si++) for (var ti = s[0][3] ? tokens.length - 1 : 0, t; t = tokens[ti]; s[0][3] ? ti-- : ti++) for (var opi = 0, op; op = s[opi]; opi++) if (t[1] == op[0]) {
                    // exception case .. ".Normalized" and ".Length" properties are re-routed (so they work on scalars etc ..)
                    if (op[2] == 2) { var arg = tokens[left()]; glue(ti - 1, ti, 5, ["Element." + op[1], "(", arg, ")"]); }
                    // unary operators (all are to the left)
                    else if (op[2]) { var arg = tokens[right()]; glue(ti, right(), 5, ["Element." + op[1], "(", arg, ")"]); }
                    // binary operators
                    else { var l = left(), r = right(), a1 = tokens[l], a2 = tokens[r]; if (op[0] == op[1]) glue(l, r, 5, [a1, op[1], a2]); else glue(l, r, 5, ["Element." + op[1], "(", a1, ",", a2, ")"]); ti--; }
                }
                return tokens;
            }
            // Glue all back together and return as bound function.
            return eval(('(' + (function f(t) { return t.map(t => t instanceof Array ? f(t) : typeof t == "string" ? t : "").join(''); })(translate(tok)) + ')'));
        }

        /**
         * 
         * @param {*} f 
         * @param {Options} options 
         */
        static graph(f, options) {
            // console.log("graph");
            return Element_graph.bind(this)(
                f,
                options,
                Element,
                tot,
                drm
            );
        }

        // webGL Graphing function. (for parametric defined objects)
        /**
         * 
         * @param {*} f 
         * @param {Options} options
         */
        static graphGL(f, options) {
            // console.log("graphGL");
            return Element_graphGL.bind(this)(
                f,
                options,
                Element,
                tot,
                counts
            );
        }

        /**
         * 
         * @param {*} f 
         * @param {Options} options 
         */
        static graphGL2(f, options) {
            // console.log("graphGL2");
            return Element_graphGL2.bind(this)(
                f,
                options,
                Element,
                tot,
                counts
            );
        }
    }

    return res;
}
