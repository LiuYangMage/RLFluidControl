import java.util.*;
import org.apache.xmlrpc.*;	
AFCCylinder test;				
SaveScalar dat;				
int Time = 50;	
int initTime = 1;
int plotTime = 48;
int picNum = 20; 
float tCurr = 0, tStep = .0075;
int simNum = 1;
String path = "D:\\RL - Dixia\\RLCYLFix\\RLCYLFix\\saved\\";
int callLearn = 16;
float Cd = 0, Cl = 0;
int resolution = 24;

XmlRpcClient client;

void setup(){
  size(512,256);//size(xLengths*resolution,yLengths*resolution); 
  try{
    client = new XmlRpcClient("http://localhost:8000");
    Vector params = new Vector();
    params.addElement(new Integer(-1));
    Object result = client.execute("init", params);
    print(result);
  } 
  catch (Exception ex){
    println(ex);
  }
  
  setUpNewSim(simNum);  
  
}

void draw() {        
  if (test.t<Time) {        
    test.update2();              
    dat.addData03(test.t, test.xi1, test.xi2, test.force, test.flow.p);  
    if (test.t>initTime){
      callLearn--;
      Cd += test.force.x;
      Cl += test.force.y;
      
      if(callLearn<=0){
        callLearn = 16;
        Cd = Cd/callLearn*2/resolution;
        Cl = Cl/callLearn*2/resolution;
        
        // call for action
        float[] xiNext = callAction(Cl, Cd);
        test.xi1 = xiNext[0];
        test.xi2 = xiNext[1];
        test.xi1_m = 5*test.xi1;
        test.xi2_m = 5*test.xi2;
      }   
      
      if(test.t>plotTime){
        picNum--;
        if(picNum <= 0){
          test.display();
          saveFrame("saved/"+str(simNum) + "/" +"frame-#######.png");
          picNum = 10;
        }
      }
    }                
  }        
  else { 
    test.display();
    saveFrame(path + "\\" + str(simNum) + "\\result.png");        
    dat.finish(); 
    
    try{
      Vector params = new Vector();
      params.addElement(new Integer(1000));
      Object result = client.execute("train", params);
      print(result);
      client.execute("save", params);
    } 
    catch (Exception ex){
      println(ex);
    }
    simNum = simNum + 1;   
    setUpNewSim(simNum);
  }        
}        
        
void keyPressed(){   // close and save everything when the space bar is pressed               
    dat.finish();        
    exit();        
}

// New Eposide
void setUpNewSim(int runNum){       
  int xLengths = 16, yLengths = 8, zoom = 100/resolution, Re = 500;//default ufree=1         
  float dR = .125, gR = .2;         
  float xi1 = 0, xi2 = 0, theta = PI/3;  
  smooth();
  
  if (zoom <= 1){zoom = 1;}
  
  test = new AFCCylinder(resolution, Re, dR, gR, theta, xi1, xi2, tStep, xLengths, yLengths,zoom, true);          
  dat = new SaveScalar("saved/"+str(runNum)+".txt", (float)resolution, (float)xLengths, (float)yLengths, 32);      
  
  new File(path + str(runNum)).mkdir();
  
  try{
    Vector params = new Vector();
    params.addElement(new Integer(-1));
    client.execute("start_episode", params);  
  }
  catch (Exception ex){
    println(ex);
  }
}

// call Action
float[] callAction(float Cl, float Cd){
  float[] XI = new float[2]; 
  Vector params = new Vector();
  params.addElement(String.valueOf(Cl) + "_" + String.valueOf(Cd));
  
  try{
    Object result = client.execute("request_stochastic_action", params);
    String[] newString = ((String)result).split("_");
    XI[0] = Float.parseFloat(newString[0]);
    XI[1] = Float.parseFloat(newString[1]);
  }
  catch (Exception ex){
    println(ex);
  }
  
  return XI;
}
