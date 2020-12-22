/**
 * 
 * @param {*} w 
 * @param {*} h 
 * @param {*} data 
 * @param {*} canvas 
 * @param {*} f
 */

export function graph_f1(w, h, data, canvas, f) {
    // one parameter function.. go over x range, use result as y.
    for (var px = 0; px < w; px++) {
        var res = f(px / w * 2 - 1);
        res = Math.round((res / 2 + 0.5) * h);
        if (res > 0 && res < h - 1)
            data.data.set(
                [0, 0, 0, 255],
                res * w * 4 + px * 4
            );
    }
}
