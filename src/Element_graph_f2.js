/**
 * 
 * @param {*} w 
 * @param {*} h 
 * @param {*} data 
 * @param {*} canvas 
 * @param {*} f
 */

export function graph_f2(w, h, data, canvas, f) {
    // two parameter functions .. evaluate for both and set resulting color.
    for (var px = 0; px < w; px++) {
        for (var py = 0; py < h; py++) {
            var res = f(px / w * 2 - 1, py / h * 2 - 1);
            res = res.buffer ? [].slice.call(res) : res.slice ? res : [res, res, res];
            data.data.set(res.map(x => x * 255).concat([255]), py * w * 4 + px * 4);
        }
    }
}
