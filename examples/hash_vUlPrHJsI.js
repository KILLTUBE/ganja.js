// a PGA implementation of Keenan Crane's Heatflow method for Geodesic distance.
// https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/
Algebra(3,0,1,()=>{ 
    
  // Parse a waverfront OBJ file into a list of vertex positions and face vertex-indices.    
  const parseOBJ = (txt,scale=10)=>{
    var t = txt.split('\n'), v=t.filter(x=>x.match("v ")),
        f = t.filter(x=>x.match("f ")).map(x=>x.split(/\s+/).slice(1,4).reverse().map(x=>parseInt(x)-1));
    return [v.map(x=>x.split(/\s+/).slice(1,4).map(x=>scale*parseFloat(x))), f]
  }
  
  // Convert a point and face list to a Half-Edge mesh with PGA entities for points, lines and planes
  const halfEdge = ([pts,faces])=>{
    var HE = []; pts = pts.map(x=>!(1e0+x*[1e1,-1e2,1e3]+1e2))
    faces = faces.map(([v1,v2,v3],i)=>{
      var face = {vtx : [v1,v2,v3], plane : pts[v1]&pts[v2]&pts[v3], idx : i};
      HE.push(...[[v1,v2],[v2,v3],[v3,v1]].map(([v1,v2],j)=>({ from:v1, to:v2, line:pts[v1]&pts[v2], face })));
      for (var j=0; j<3; ++j) face.edge = HE[i*3+((j+2)%3)].next = HE[i*3+((j+1)%3)].prev = HE[i*3+j];
      return face;
    });
    pts = pts.map((x,i)=>({idx:i,vtx:x,edge:HE.find(x=>x.from==i)}))
    HE.forEach(E=>E.twin = HE.find(E2=>(E2.to==E.from)&&(E.to==E2.from)));
    const Edges = (p,f,s)=>{ var e=p.edge; while(1){ s=f(e,s); e=(p.plane)?e.next:e.twin.next; if (e==p.edge) return s; }}
    return [pts,HE,faces,Edges];
  }
  
  // Implement the heatFlow method.
  const heatFlow = ([pts,HE,faces,Edges],vtx=[198])=>{
    // The cotan formula in PGA directly uses PGA lines and planes.
    var cotan = h=>(h.next.line|h.prev.line).s/h.face.plane.Length;
    
    // Build the Mass matrix, a diagonal matrix with barycentric dual weight per vertex.
    var M = [...Array(pts.length)].map(x=>[...Array(pts.length)].map(x=>0)); 
    pts.forEach((p,i)=>Edges(p,e=>{ M[i][i] += e.face.plane.Length/6 }));
    
    // The Laplacian matrix in its positive-definite version.
    var L = [...Array(pts.length)].map((x,i)=>[...Array(pts.length)].map((x,j)=>i==j?1E-8:0));
    pts.forEach((p,i)=>Edges(p,e=>{ L[i][i] -= L[i][e.to] = -(cotan(e) + cotan(e.twin))/2; }));
    
    // The flow matrix, mass + laplacian*average_edge_length
    var F = M + L*(HE.reduce((s,e)=>s+e.line.Length,0)/HE.length)**2;
    
    // Our vector with initial heat per vertex, Integrated one timestep.
    var heat = [...Array(pts.length)].map(x=>0); vtx.forEach(v=>{heat[v]=1});
    let u = heat/F;              

    // The gradient (normalized) and divergence of the heat flow u.
    var X = faces.map(f=>Edges(f,(e,grad)=>grad+f.plane.Normalized|e.line*u[e.prev.from],0).Normalized);
    var div = pts.map((p,i)=>0.5*Edges(p,(h,sum)=>
          sum + cotan(h)*(h.line^X[h.face.idx]).e123 + cotan(h.prev)*(h.prev.twin.line^X[h.face.idx]).e123
        ,0));  
        
    // Now solve the poisson equation to find the distance.    
    let phi = div/L, min = Math.min(...phi), max = Math.max(...phi);
    return [pts,faces,phi.map(x=>(x-min)/(max-min))];
  }

  // Fetch, process and display the bunny.
  fetch("https://enki.ws/bunny_low.obj").then(x=>x.text())
  .then(parseOBJ)
  .then(halfEdge)
  .then(heatFlow)
  .then(([pts,faces,phi])=>{
    var obj = {data: pts.map(x=>x.vtx), idx: faces.map(f=>f.vtx).flat(), color: phi.map(x=>[1-x*1.4,Math.max(0,1-x*4),Math.max(0,1-x*8)]).flat()};
    document.body.appendChild(this.graph([obj],{camera:1-1e03,width:'100%',height:'100%',gl:1,
    shader:`gl_FragColor =  gl_FragColor*0.2 + 0.4*l*vec4(Col,1.) + 0.7 * vec4(Col,1.);
            float m = mod(Col.r*40.,1.); 
            gl_FragColor += m<0.8?0.:0.3*sin((m-0.8)*5.*3.1415);`}));
  })
});