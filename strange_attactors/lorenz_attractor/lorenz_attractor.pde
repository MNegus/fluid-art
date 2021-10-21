// Daniel Shiffman
// http://codingtra.in
// http://patreon.com/codingtrain
// Code for this video: https://youtu.be/f0lkz2gSsIk

import peasy.*;
import processing.svg.*;

boolean record;

float x = 0.01;
float y = 0;
float z = 0;

float a = 20;
float b = 28;
float c = 8.0/3.0;

float dt = 0.005;
float tmax = 40;
float t = 0;

float x_ang = 0;
float y_ang = 0;

ArrayList<PVector> points = new ArrayList<PVector>();

//PeasyCam cam;

void setup() {
  size(800, 600, P2D);
  colorMode(HSB);
  //cam = new PeasyCam(this, 500);
  noFill();
  // Create points
  while (t <= tmax) {
    
    float dx = (a * (y - x))*dt;
    float dy = (x * (b - z) - y)*dt;
    float dz = (x * y - c * z)*dt;
    x = x + dx;
    y = y + dy;
    z = z + dz;
  
    points.add(new PVector(x, y, z));
    
    t = t + dt;
  }
}

void draw() {
  if (record) {
    //beginRaw(SVG, "output.svg");
    String filename = "a_" + str(a) + "-b_" + str(b) + "-c_" + str(c) + ".svg";
    noFill();
    beginRecord(SVG, filename);
  }
  
  //background(255);
  
  //translate(0, 0, -80);
  translate(width/2, height/2);
  scale(5);
  stroke(0);
  strokeWeight(0.25);
  noFill();

  beginShape();
  for (PVector v : points) {
    //vertex(v.x, v.y,v.z);
    
    // Loads in (x, y)
    float x = v.x; 
    float y = v.y;
    float z = v.z;
    
    // Rotates (x, y)
    float rx = cos(y_ang) * x + sin(y_ang) * y;
    float ry = -sin(x_ang) * sin(y_ang) * x + cos(x_ang) * y - sin(x_ang) * cos(y_ang) * z;
    
    //float rx = x;
    //float ry = y;
    
    vertex(rx, ry);
    //vertex(x, y, z);
    //PVector offset = PVector.random3D();
    //offset.mult(0.1);
    //v.add(offset);
  }
  endShape();
  
  //noLoop();
  

  if (record) {
    //endRaw();
    endRecord();
    record = false;
  }

  //println(x,y,z);
}

// Hit 'r' to record a single frame
void keyPressed() {
  if (key == ' ') {
    //clear();
    x_ang = random(1) * PI;
    y_ang = random(1) * PI;
  }
  if (key == 'r') {
    record = true;
  }
}
