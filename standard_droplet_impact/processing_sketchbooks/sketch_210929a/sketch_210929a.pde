import processing.svg.*;

int timestep;
int inc;
int maxtimestep;
int frameno;
float xmax;
float ymax;
boolean record;

void setup() {
  //size(512, 512);
  size(180, 180);
  background(255);
  
  xmax = 2;
  ymax = 2 * xmax;
  frameno = 1;
  timestep = 0;
  inc = 37;
  maxtimestep = 1500;
  record = true;
  
}


void draw() {
  clear();
  background(255);
  noFill();
  stroke(0);
  
  frameRate(20);
  
  if (record) {
    String filename = "frames/frame_" + frameno + ".svg";
    noFill();
    beginRecord(SVG, filename);
  }
  
  // Draws border
  noFill();
  square(0, 0, width);
  
  // Draws interface
  String data_name = "interface_data/interface_" + str(timestep) + ".txt";
  String[] lines = loadStrings(data_name);
  noFill();
  beginShape(); 
  for (int i = 0 ; i < lines.length; i++) {
    String[] pieces = split(lines[i], ',');
    float x = float(pieces[0]) ;
    float y = float(pieces[1]);
    if ((x <= xmax) && (x >= -xmax) && (y <= ymax)) {
      float px = (x+ xmax) * width / (2 * xmax);
      float py = 0.8 * height - y * height / ymax;
      vertex(px, py);
    }
  }
  endShape();
  
  
  
  // Saves
  if (record) {
    endRecord();
    frameno++;
  }
  
  timestep += inc;
  
  if (timestep >= maxtimestep) {
    noLoop();
  }
}

//String[] lines = loadStrings("interface_200.txt");
//println("there are " + lines.length + " lines");
