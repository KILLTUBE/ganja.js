https://enkimute.github.io/ganja.js/examples/coffeeshop.html#NpsKhuwlO
// Ganja.js does not handle R2 or R3 natively, but this is easy to work around.
// Lets do it for R2.

Algebra(2,()=>{

  // R2 is often looked at point-based. That means vectors are points...    
  // In this space we can not distinguish vectors and points.. 

  var point  = (x,y)=> x*1e1 + y*1e2;
  var vector = point;        
  
  // It only has rotations around the origin.
  var rotate = (point, alpha)=> Math.E**( alpha/2 * 1e12 ) >>> point;
  
  // Translations are handled by adding vectors
  var translate = (point, translation_vector)=> point + translation_vector;
  
  // Create some points
  var [A,B,C] = [point(0,0), point(1,0.5), point(0.4,-0.4)];
  
  // Rotate then translate them
  var [D,E,F] = translate(rotate([A,B,C], Math.PI/4 ), vector(0.5,-0.5));

  // Now to graph them we quickly go to 2D PGA to graph.
  // Just ignore this block basically.
  var graph = Algebra(2,0,1).inline((ar,opts)=>{
      var upgrade_to_pga = x=>{
        while (x.call) x = x.call();
        if (x instanceof Array) return x.map(upgrade_to_pga);
        if (x instanceof Float32Array) return new Element([0,0,0,0,-x[2],-x[1],1,0])
        return x;
      }
      return this.graph(upgrade_to_pga(ar),opts)
  })    
    
  // Now graph these points, polygons, lines between them etc.
  document.body.appendChild(graph([
    0xFF0000, A,"A",B,C, [A,B,C],
    0x0000FF, D,"D",E,F, [D,E],[E,F],[D,F]
  ]))    
 
    
    
    
})