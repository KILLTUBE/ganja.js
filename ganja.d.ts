interface Algebra {

}

declare class MultiVector {
    constructor(a: any);
    Grade(grade: any, res: any): any;
}

interface Element {

}

interface Options {
    scale?: number;
    conformal?: boolean;
    even?: boolean;
    dual?: number;
    camera?: Element;
    width?: number;
    height?: number;
    canvas?: HTMLCanvasElement;
    ni?: Element;
    no?: Element;
    gl?: boolean; // WebGL
    up?: boolean; // WebGL, return graphGL2
    alpha?: boolean; // WebGL
    spin?: number; // WebGL camera
    h?: number; // WebGL camera
    p?: number; // WebGL camera
    z?: number; // WebGL camera
    posx?: number; // WebGL camera
    posy?: number; // WebGL camera
    posz?: number; // WebGL camera
    maxSteps?: number; // WebGL, default 64
    stepSize?: number; // WebGL, default 0.25
    thresh?: number; // WebGL, default 0.2
    devicePixelRatio?: number; // default 1
    pointRadius?: number; // default 1
    lineWidth?: number; // default 1, doesn't work anyway (WebGL limitation)
    /*
        shader: `
            gl_FragColor =  gl_FragColor*0.2 + 0.4*l*vec4(Col,1.) + 0.7 * vec4(Col,1.);
            float m = mod(Col.r*40.,1.); 
            gl_FragColor += m<0.8?0.:0.3*sin((m-0.8)*5.*3.1415);
        `
    */
    shader?: string; // default '', special shader for programcol
    ipns?: boolean;
    still?: boolean;
    animate?: boolean;
    useUnnaturalLineDisplayForPointPairs?: boolean;
    htmlText?: boolean;
    grid?: boolean;
    metric?: number[];
    Cayley?: string[][];
    basis?: string[];
    q?: number;
    r?: number;
    mix?: boolean;
    graded?: boolean;
    grades?: number[];
    baseType?: Array<number> | Float32Array | Float64Array;
    tot?: number;
    over?: Algebra;
    fontSize: number;
    arrows: boolean; // graph with arrows or not (looks a bit convoluted in some examples)
    format?: "table" | "canvas";
}

type ElementInput  = Element             | Array<any> | string | number;
type ElementOutput = Element | Element[] | Array<any> | string | number;
interface GalculatorButtonObject {
    color?: string;
    label: string;
    click: ()=>void;
    help?: string;
    el?: HTMLDivElement;
}

// Fix TS
declare function parseFloat(x: number | string): number;

interface Number {
    call: any;
}

//import * as D3 from "d3";
//declare function d3(string: any): D3;
