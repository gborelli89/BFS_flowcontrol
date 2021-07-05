#include <Stepper.h> 
 
const int stepsPerRevolution = 500;
String v1;
String v2; 
  
//Pins
Stepper motorY1(stepsPerRevolution, 27,31,29,33); 
Stepper motorY2(stepsPerRevolution, 26,30,28,32);
Stepper motorY3(stepsPerRevolution, 35,39,37,41);
Stepper motorX1(stepsPerRevolution, 34,38,36,40);
Stepper motorX2(stepsPerRevolution, 43,47,45,49);
Stepper motorX3(stepsPerRevolution, 42,46,44,48); 

void setup() 
{   
    Serial.begin(9600);
    //Determina a velocidade inicial do motor 
    motorX1.setSpeed(20);
    motorY1.setSpeed(20);
    motorX2.setSpeed(20);
    motorY2.setSpeed(20);
    motorX3.setSpeed(20);
    motorY3.setSpeed(20);
} 
  
void loop() 
{ 
    String motorID = "0";

    if(Serial.available()){
        
        String v1 = Serial.readStringUntil(';'); 
        String v2 = Serial.readStringUntil('\0');
        
        motorID = v1;
        int motorSteps = v2.toInt();

        if(motorID=="x1"){
          motorX1.step(motorSteps);
        }
        if(motorID=="y1"){
          motorY1.step(motorSteps);
        }
        if(motorID=="x2"){
          motorX2.step(motorSteps);
        }
        if(motorID=="y2"){
          motorY2.step(motorSteps);
        }
        if(motorID=="x3"){
          motorX3.step(motorSteps);
        }
        if(motorID=="y3"){
          motorY3.step(motorSteps);
        }
 
     }

     Serial.println(motorID);
     delay(500);
}
