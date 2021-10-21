import peasy.*;

// Initial point
float x = 0.01;
float y = 0.01;
float z = 0.01;

// Time parameters
float dt = 0.0005;
float tmax = 500;
float t = 0;

// Coefficients
float[] coeffs = new float[27];

ArrayList<PVector> points = new ArrayList<PVector>();

PeasyCam cam;

void setup() {
  // Canvas and camera set up
  size(800, 600, P3D);
  colorMode(HSB);
  cam = new PeasyCam(this, 500);
  noFill();
  
  // Initialise all coefficients to 0
  for (int i = 0; i < 27; i = i + 1) {
    coeffs[i] = 0;
  }
  
  // Lorenz attractor
  //float sigma = 10;
  //float beta = 8.0 / 3.0;
  //float rho = 28;
  //coeffs[0] = - sigma;
  //coeffs[1] = sigma;
  //coeffs[9] = rho;
  //coeffs[10] = -1;
  //coeffs[13] = -1;
  //coeffs[20] = - beta;
  //coeffs[21] = 1;
  
  // Four-Wing Attractor
  float a = 2;
  float b = -0.1;
  float c = 10;
  float d = -4;
  float e = -10;
  float f = -10;
  coeffs[0] = a;
  coeffs[5] = c;
  coeffs[9] = b;
  coeffs[10] = d;
  coeffs[13] = -1;
  coeffs[20] = e;
  coeffs[21] = f;
  
  
  // Create points
  while (t <= tmax) {
    
    float dx = (coeffs[0] * x + coeffs[1] * y + coeffs[2] * z + coeffs[3] * x * y + coeffs[4] * x * z + coeffs[5] * y * z + coeffs[6] * x * x + coeffs[7] * y * y + coeffs[8] * z * z)*dt;
    float dy = (coeffs[9] * x + coeffs[10] * y + coeffs[11] * z + coeffs[12] * x * y + coeffs[13] * x * z + coeffs[14] * y * z + coeffs[15] * x * x + coeffs[16] * y * y + coeffs[17] * z * z)*dt;
    float dz = (coeffs[18] * x + coeffs[19] * y + coeffs[20] * z + coeffs[21] * x * y + coeffs[22] * x * z + coeffs[23] * y * z + coeffs[24] * x * x + coeffs[25] * y * y + coeffs[26] * z * z)*dt;
    x = x + dx;
    y = y + dy;
    z = z + dz;
  
    points.add(new PVector(x, y, z));
    
    t = t + dt;
    println(t);
    //println(x, y, z);
  }
}

void draw() {
  
  background(255);
  
  translate(0, 0, -80);
  //translate(width/2, height/2);
  scale(20);
  stroke(0);
  strokeWeight(0.1);
  noFill();

  beginShape();
  for (PVector v : points) {
    
    // Loads in (x, y)
    float x = v.x; 
    float y = v.y;
    float z = v.z;
    
    vertex(x, y, z);
  }
  endShape();
  
  //noLoop();
  
}
