/**********************************
 SaveScalar class
 
 Saves scalar data to a text file with customizable header
 
 example code:
 // dat = new SaveData("pressure.txt",test.body.);
 dat.addData(test.t, test.flow.p);
 dat.finish();
 ***********************************/

class SaveScalar{
  PVector force;
  PrintWriter output;
  int n;
  float  resolution, xi1, xi2, theta, D;
  float centX, centY;
  int numTheta;
  
  SaveScalar(String name){
    output = createWriter(name);
    output.println("%% Force coefficients using processing viscous simulation");
    output.println();
    output.println("%% Fellowing: t, force.x, force.y");
    output.println();
  }
  
  SaveScalar(String name, float res, float xLen, float yLen, int num){
    output = createWriter(name);
    output.println("%% Force and pressure coefficients using processing viscous simulation");
    output.println();
    output.println("%% Fellowing: t, force.x, force.y");
    output.println();
    this.resolution = res;
    this.D = res;
    float n = xLen*res;
    float m = yLen*res;
    this.numTheta = num;
    this.centX = n/4;
    this.centY = m/2;
  }
     
  void addData(float t,PVector force){
    output.println(t + " " + force.x + " " +force.y);
    output.println("");
  }
  
  void addData02(float t,PVector force, Field pres){
    output.print(t + " " + force.x + " " +force.y + " ");
    for(int i=0;i<numTheta; i++){
      float xPre = cos((float)i/numTheta*PI*2)*D/2 + centX;
      float yPre = sin((float)i/numTheta*PI*2)*D/2 + centY;
      
      float pdl = pres.linear( xPre, yPre );
      output.print(pdl + " ");
    }
    //output.println("");
    output.println("");
  }
  
  void addData03(float t, float r1, float r2, PVector force, Field pres){
    output.print(t + " " + force.x + " " +force.y + " " + r1 + " " + r2 + " ");
    for(int i=0;i<numTheta; i++){
      float xPre = cos((float)i/numTheta*PI*2)*D/2 + centX;
      float yPre = sin((float)i/numTheta*PI*2)*D/2 + centY;
      
      float pdl = pres.linear( xPre, yPre );
      output.print(pdl + " ");
    }
    //output.println("");
    output.println("");
  }
  
  
  void finish(){
    output.flush(); // Writes the remaining data to the file
    output.close(); // Finishes the file
  }
} 
  
