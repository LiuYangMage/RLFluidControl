class AFCCylinder {
  BDIM flow;
  BodyUnion body;
  boolean QUICK = true, order2 = true;
  int n, m, out, up, resolution,NT=1;
  float dt, t, D, xi1, xi2,theta;//initial dt=1 only determine BDIM's update u.
  float xi1_m, xi2_m, theta_m, dR, gR, r, dphi1,dphi2;
  FloodPlot flood;
  PVector force;

  AFCCylinder (int resolution, int Re,float dR, float gR,  float theta,  float xi1,  float xi2, float dtReal, int xLengths, int yLengths,  int zoom, boolean isResume) {
    n = xLengths*resolution;
    m = yLengths*resolution;
    this.resolution = resolution;
    this.xi1 = xi1;
    this.xi2 = xi2;
    this.dR = dR;
    this.gR = gR;
    this.theta = theta;
    this.dt = dtReal*this.resolution;
    xi1_m=5*xi1;
    xi2_m=5*xi2;
    theta_m=theta;


    Window view = new Window(0, 0, n, m); // zoom the display around the body
    D=resolution;

    float r=(D/2+gR*D+dR*D/2);
    body =new BodyUnion(new CircleBody(n/4, m/2, D, view),
    new CircleBody(n/4+r*cos(theta_m), m/2-r*sin(theta_m), dR*D, view),
    new CircleBody(n/4+r*cos(theta_m), m/2+r*sin(theta_m), dR*D, view));
    
    flow = new BDIM(n,m,dt,body,(float)D/Re,QUICK);
    if(isResume){
      flow.resume("saved\\init\\init.bdim");
    }
    
    flood = new FloodPlot(view);
    flood.range = new Scale(-1, 1);
    flood.setLegend("vorticity"); 
  }


 void update2(){
   //dt = flow.checkCFL();
   flow.dt = dt;
   dphi1 = (2*xi1_m*dt)/(dR*D);
   dphi2 = (2*xi2_m*dt)/(dR*D);
   body.bodyList.get(1).rotate(dphi1);//change index try
   body.bodyList.get(2).rotate(dphi2);
              
   flow.update(body);
   if (order2) {flow.update2(body);}
   print("t="+nfs(t,2,2)+";  ");
   t += dt/resolution;  //nonedimension
  
   force = body.bodyList.get(0).pressForce(flow.p).mult(-1);
   print("drag="+nfs(force.x*2/D, 2, 2)+";  ");
   println("lift="+nfs(force.y*2/D, 2, 2)+";  ");
 }

 void update(){
    for ( int i=0 ; i<NT ; i++ ) {
      if (flow.QUICK) {
        dt = flow.checkCFL();
        flow.dt = dt;
      }
       
       dphi1 = (2*xi1_m*dt)/(dR*D);
       dphi2 = (2*xi2_m*dt)/(dR*D);
       body.bodyList.get(1).rotate(dphi1);//change index try
       body.bodyList.get(2).rotate(dphi2);
              
       flow.update(body);
       if (order2) {flow.update2(body);}
       print("t="+nfs(t,2,2)+";  ");
       t += dt/resolution;  //nonedimension
  
      force = body.bodyList.get(0).pressForce(flow.p).mult(-1);
      print("drag="+nfs(force.x*2/D, 2, 2)+";  ");
      println("lift="+nfs(force.y*2/D, 2, 2)+";  ");   
  }
}

  void display() {
    flood.display(flow.u.curl());
    body.display();
    flood.displayTime(t);
  }
}
