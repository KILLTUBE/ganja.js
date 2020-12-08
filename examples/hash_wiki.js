// https://discourse.bivector.net/t/eigenblades-of-linear-functions/276

Algebra(3,0,1,()=>{
    var {E,PI}=Math;
    var BV = 1e12 + .3e03;
    var d = {motor:E**BV, xRange:2, yRange:2, zRange:2};
    var orbit_1d = s => E**((s-0.5)*5*BV) * (1+.1e01)
    document.body.appendChild(this.graph(()=>{
      return [
          d,
          0x0000ff,orbit_1d, 
          0x008800,BV 
      ];
    },{animate:1,grid:0,gl:1,h:-1,camera:0e1,alpha:1})).style.backgroud='transparent';
});

disp=document.querySelector("svg")||document.querySelector("canvas");
if (disp && disp.value) {
    disp.style.width=disp.style.height="100%";
    if (disp.update && (!disp.options || !disp.options.animate)) {
        disp.update(disp.value);
    }
} 