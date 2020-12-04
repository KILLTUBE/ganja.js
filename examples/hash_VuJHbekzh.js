// https://enkimute.github.io/ganja.js/examples/coffeeshop.html#VuJHbekzh
// Create a Clifford Algebra with 2,0,1 metric. 
Algebra(2,0,1,()=>{

  var {PI,E,abs} = Math;

 
  var arrow = (a,b)=>{
    var o = !1e0;
    var line = a&b;
    var xaxis = 1e1;
    var len = line.Length;
    var l2 = 0.1;
    var arrow_x = [[o,o + !(len-l2)*1e1],
                   [o + !(l2*.5e2+(len-l2)*1e1),o + !len*1e1,o - !(l2*.5e2-(len-l2)*1e1)]];
    var tran = (1+a/o).Normalized;
    var l1 = ~tran >>> line.Normalized;
    var rot = (1 - l1*1e2).Normalized;
    return (tran*rot)>>>arrow_x;  
  }
    
    var p1 = !(1e0 - 1e1);
    var p2 = !(1e0 + 1e1 + 0.925e2);
  var canvas = document.body.appendChild(this.graph(()=>{

    var orig = !(1e0-1.5e1 + 1.5e2);
    var l = (1e1+.5e2).Normalized;
    
    var p12 = (p1.Normalized)-orig+p2.Normalized;
    return [
      0xFF0000,p1,"a",...arrow(orig,p1),
      0x0000FF,p2,"b",...arrow(orig,p2),
      0x999900,p12,"a+b",0xCCCCFF,...arrow(p1,p12), 0xFFCCCC,...arrow(p2,p12),
      0xFFAAFF,(p2|(orig&p1))*(orig&p1),"((b&bullet;a)/|a|²)a",[(p2|(orig&p1))*(orig&p1),p2],...arrow(orig,((p2|(orig&p1))*(orig&p1)).Normalized),
      0x00FF00,l,0x008800,"line",...arrow(orig,(orig|l)*l ),...arrow((orig|l)*l, ((orig|l)*l + !1e0)/3 ),
      0,1e1+1.5e0,1e2-1.5e0
    ];  
  },{grid:true, lineWidth:4, fontSize:1.5, pointRadius:1, animate:true, alpha:true,gl:0}));
  
  canvas.style.background='transparent';

});