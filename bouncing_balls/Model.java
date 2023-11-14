package bouncing_balls;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.lang.Math.*;
import java.util.Vector;


class Model {

	double areaWidth, areaHeight;
	double PI = 3.14159265;
	
	Ball [] balls;
	Rectangle [] rectangles;
	HashMap<String, Boolean> collisionEligible = new HashMap<>(); //This hashMap applies only for collisions between Rectangles and Balls

	Model(double width, double height) {
		areaWidth = width;
		areaHeight = height;
		
		// Initialize the model with a few balls.
		balls = new Ball[4];
		rectangles = new Rectangle[1];

		balls[0] = new Ball(width / 3, height * 0.9, 1.2, 1.6, 0.2, 1);
		balls[1] = new Ball(2 * width / 3, height * 0.7, -0.6, 0.6, 0.282, 2);
		balls[2] = new Ball(0.7 * width / 3, height * 0.1, -0.3, 0.3, 0.10, 0.5);
		balls[3] = new Ball(1.4 *width / 3, height * 0.8, 1.1, 1.5, 0.2, 1);
		//balls[4] = new Ball(1.5 * width / 3, height * 0.2, -1.6, 3.6, 0.682, 3);

		rectangles[0] = new Rectangle(width / 2.5,height * 0.4,width / 2.5,height * 0.4,(width / 2.5)+0.25, height * 0.4,width / 2.5,height * 0.4+0.25,
				(width / 2.5)+0.25,height * 0.4+0.25,1.2,1.6,PI/4,1, PI/4,0.25,0.25);

		//enables all rectangles to collide with all balls
		for (int i = 0; i<rectangles.length; i++){
			for (int j = 0; j<balls.length; j++) {

				int touple[] = new int[2];
				touple[0] = i;
				touple[1] = j;

				collisionEligible.put(Arrays.toString(touple), true);
			}
		}
	}

	void step(double deltaT) {		

		// Saves the old information about the balls.

		// Checks if the balls will overlap each other in the next frame and if so simulate a collision.
		checkBalls(deltaT);

		checkRectangles(deltaT);
		for (int j = 0; j < balls.length; j++){


			// Detect collision with the border.
			if (balls[j].x < balls[j].radius || balls[j].x > areaWidth - balls[j].radius) {
				balls[j].vx *= -1; // Change direction of ball.
				for (int i = 0; i < rectangles.length; i++){
					int touple[] = new int[2];
					touple[0] = i;
					touple[1] = j;
					collisionEligible.replace(Arrays.toString(touple),true);
				}
			}
			if (balls[j].y < balls[j].radius || balls[j].y > areaHeight - balls[j].radius) {
				balls[j].vy *= -1;
				for (int i = 0; i < rectangles.length; i++){
					int touple[] = new int[2];
					touple[0] = i;
					touple[1] = j;
					collisionEligible.replace(Arrays.toString(touple),true);
				}
			}
			// Compute new position according to the speed of the ball.
			double oldVy = balls[j].vy;
				balls[j].vy -= deltaT * 9.82;

				balls[j].x += deltaT * balls[j].vx;
				balls[j].y += deltaT * (balls[j].vy + oldVy) / 2; // Computes the new positions with the average speed over the timeframe.

		}

		for (int i = 0; i < rectangles.length; i++){


				//Detect collision with the border.
				if (rectangles[i].x1 <= 0 || rectangles[i].x2 <= 0 || rectangles[i].x3 <= 0 || rectangles[i].x4 <= 0 ||
						rectangles[i].x1 >= areaWidth || rectangles[i].x2 >= areaWidth || rectangles[i].x3 >= areaWidth ||
					rectangles[i].x4 >= areaWidth) {
					rectangles[i].vx *= -1;
					rectangles[i].rotation *= -1;
					for (int j = 0; j < balls.length; j++){
						int touple[] = new int[2];
						touple[0] = i;
						touple[1] = j;
						collisionEligible.replace(Arrays.toString(touple),true);
					}
					//Glöm inte att fixa rotationen
				}

				if (rectangles[i].y1 <= 0 || rectangles[i].y2 <= 0 || rectangles[i].y3 <= 0 || rectangles[i].y4 <= 0 ||
					rectangles[i].y1 >= areaHeight || rectangles[i].y2 >= areaHeight || rectangles[i].y3 >= areaHeight ||
					rectangles[i].y4 >= areaHeight) {
					rectangles[i].vy *= -1;
					rectangles[i].rotation *= -1;
					for (int j = 0; j < balls.length; j++){
						int touple[] = new int[2];
						touple[0] = i;
						touple[1] = j;
						collisionEligible.replace(Arrays.toString(touple),true);
					}
					//Glöm inte att fixa rotationen
				}

				// Compute new position according to the speed of the rectangle.
				double oldVy = rectangles[i].vy;
				rectangles[i].vy -= deltaT * 9.82;

				//update the four corners.
				rectangles[i].x1 += deltaT * rectangles[i].vx;
				rectangles[i].y1 += deltaT * (rectangles[i].vy + oldVy) / 2;
				rectangles[i].x2 += deltaT * rectangles[i].vx;
				rectangles[i].y2 += deltaT * (rectangles[i].vy + oldVy) / 2;
				rectangles[i].x3 += deltaT * rectangles[i].vx;
				rectangles[i].y3 += deltaT * (rectangles[i].vy + oldVy) / 2;
				rectangles[i].x4 += deltaT * rectangles[i].vx;
				rectangles[i].y4 += deltaT * (rectangles[i].vy + oldVy) / 2;

				//taking the rotation into consideration
				double centerX = (Math.abs(rectangles[i].x1 + rectangles[i].x4))/2;
				double centerY = (Math.abs(rectangles[i].y1 + rectangles[i].y4))/2;

				//moves the rectangle to the origin
				double x1o = rectangles[i].x1 - centerX;
				double y1o = rectangles[i].y1 - centerY;

				double x2o = rectangles[i].x2 - centerX;
				double y2o = rectangles[i].y2 - centerY;

				double x3o = rectangles[i].x3 - centerX;
				double y3o = rectangles[i].y3 - centerY;

				double x4o = rectangles[i].x4 - centerX;
				double y4o = rectangles[i].y4 - centerY;

				//rotate the rectangle
				double deltaAngle = rectangles[i].rotation * deltaT;
				double rotatedX1 = x1o * Math.cos(deltaAngle) - y1o * Math.sin(deltaAngle);
				double rotatedY1 = x1o * Math.sin(deltaAngle) + y1o * Math.cos(deltaAngle);

				double rotatedX2 = x2o * Math.cos(deltaAngle) - y2o * Math.sin(deltaAngle);
				double rotatedY2 = x2o * Math.sin(deltaAngle) + y2o * Math.cos(deltaAngle);

				double rotatedX3 = x3o * Math.cos(deltaAngle) - y3o * Math.sin(deltaAngle);
				double rotatedY3 = x3o * Math.sin(deltaAngle) + y3o * Math.cos(deltaAngle);

				double rotatedX4 = x4o * Math.cos(deltaAngle) - y4o * Math.sin(deltaAngle);
				double rotatedY4 = x4o * Math.sin(deltaAngle) + y4o * Math.cos(deltaAngle);

				//moves the rectangle back to it's original position
				rectangles[i].x1 = rotatedX1 + centerX;
				rectangles[i].y1 = rotatedY1 + centerY;
				rectangles[i].x2 = rotatedX2 + centerX;
				rectangles[i].y2 = rotatedY2 + centerY;
				rectangles[i].x3 = rotatedX3 + centerX;
				rectangles[i].y3 = rotatedY3 + centerY;
				rectangles[i].x4 = rotatedX4 + centerX;
				rectangles[i].y4 = rotatedY4 + centerY;
		}
		
	}

	private void checkRectangles(double deltaT) {

		HashMap< String, Coordinate> collisionPointMap = getRectangleCollision();
		//calculate the distance
		for (int i = 0;  i<rectangles.length; i++) {
			for (int j = 0; j < balls.length; j++) {
				int touple[] = new int[2];
				touple[0] = i;
				touple[1] = j;
				//Copy old pair of rectangle and ball
				Rectangle oldRectangle = new Rectangle(rectangles[i]);
				Ball oldBall = new Ball(balls[j]);

				double rectangleDistance = Math.sqrt(Math.pow
						((collisionPointMap.get(Arrays.toString(touple)).getX()),2) +
						Math.pow(collisionPointMap.get(Arrays.toString(touple)).getY(),2)
				);

				boolean bol = collisionEligible.get(Arrays.toString(touple));
				if (rectangleDistance <= balls[j].radius){
					if(bol){

					double centerX =  rectangles[i].x1 + (rectangles[i].x4 - rectangles[i].x1)/2;
					double centerY = rectangles[i].y1 + (rectangles[i].y4 - rectangles[i].y1)/2;

					double collisionNormalX = centerX - balls[j].x;
					double collisionNormalY = centerY - balls[j].y;

					double collisionNormalVectorLength = Math.sqrt(Math.pow(collisionNormalX, 2) + Math.pow(collisionNormalY, 2));

					//Normalizes the vector
					collisionNormalX /= collisionNormalVectorLength;
					collisionNormalY /= collisionNormalVectorLength;

					//calculate the Linear Impulse

					//u1 initial speed of rectangle along the collision normal
					//u2 initial speed of ball along the collision normal
					//v1 speed of rectangle along the collision normal after collision
					//v2 speed of ball along the collision normal after collision
					//m1 is the mass for the rectangle, the same applies for ball with m2
					//mom is the momentum
					//rv is the relative velocity before the collision

					//Computes the dot product (projection of the speed onto the line between the centerpoints).

					double u1 = (collisionNormalX * rectangles[i].vx) + (collisionNormalY * rectangles[i].vy);
					double u2 = (collisionNormalX * balls[j].vx) + (collisionNormalY * balls[j].vy);

					double m1 = rectangles[i].mass;
					double m2 = balls[j].mass;

					double mom = m1 * u1 + m2 * u2;
					double rv = u2 - u1;

					//these formulas are for the new velocities after the collision while conserving the momentum and the kinetic energy.
					double v1 = (mom + m2 * rv) / (m1 + m2);
					double v2 = v1 - rv;

					//Calculate the Torque
					double rx = centerX - collisionPointMap.get(Arrays.toString(touple)).getX();
					double ry = centerY - collisionPointMap.get(Arrays.toString(touple)).getY();

					double e = 1.0; // coefficient of restitution for a perfectly elastic collision
					double dotProduct = (collisionNormalX * (rectangles[i].vx - balls[j].vx) + (collisionNormalY * (rectangles[i].vy - balls[j].vy)));

					double impulseMagnitude = -(1 + e) * dotProduct / (1/m2 + 1/m1);
					double impulseX = impulseMagnitude * collisionNormalX;
					double impulseY = impulseMagnitude * collisionNormalY;
					double torque = rx * impulseY - ry * impulseX;

					//Calculate the Angular Impulse
					double momentOfInertia = (1.0/12.0) * rectangles[i].mass * (rectangles[i].width * rectangles[i].width + rectangles[i].height * rectangles[i].height);
					double deltaRotation = torque * deltaT / momentOfInertia;

					//update the angular Velocity (rotation) of the Rectangle
					rectangles[i].rotation += deltaRotation;

					rectangles[i].vx += collisionNormalX * (v1 - u1);
					rectangles[i].vy += collisionNormalY * (v1 - u1);
					balls[j].vx += collisionNormalX * (v2 - u2);
					balls[j].vy += collisionNormalY * (v2 - u2);

					// This method call will provide us with the distance between the balls if they continue to move with the same speed and direction for one more frame.
					double newDistance = nextRectBallDistance(deltaT, centerX, centerY, i, j);
					double currentDistance = currentRectBallDistance(deltaT, centerX, centerY, i, j);
					// This if-statement will hinder the balls from colliding a second time when they are moving away from each other after the collision but are still
					// overlapping. If the balls collide a second time, which sould not happend, then the balls will bounce towards each other.
					if (newDistance < currentDistance) {
						// Switching back to the old speed and direction makes sure they will continue to move away from each other after the collision.
						rectangles[i] = oldRectangle;
						balls[j] = oldBall;
					}
					for (int k = 0;  k<rectangles.length; k++) {
						int touple2[] = new int[2];
						touple2[0] = k;
						touple2[1] = j;
						collisionEligible.replace(Arrays.toString(touple2),true);
					}

					for (int l = 0; l < balls.length; l++) {
						int touple3[] = new int[2];
						touple3[0] = i;
						touple3[1] = l;
						collisionEligible.replace(Arrays.toString(touple3),true);
					}
					collisionEligible.replace(Arrays.toString(touple),false);
				}
				}
			}
		}
	}

	private double currentRectBallDistance(double deltaT, double centerX, double centerY, int i, int j) {
		double rectX = centerX;
		double rectY = centerY;
		double ballX = balls[j].x;
		double ballY = balls[j].y;
		return Math.sqrt(Math.pow((rectX - ballX),2) + Math.pow((rectY - ballY),2));
	}

	//This method calculates the shortest distance between the rectangle and each ball
	private HashMap<String, Coordinate> getRectangleCollision() {

		HashMap<String, Coordinate> collisionPoints = new HashMap<>();

		for (int i = 0; i<rectangles.length; i++){
			for (int j = 0; j<balls.length; j++){

				//the identityvectors for the new basis
				double iV1x = Math.cos(rectangles[i].angle);
				double iV1y = Math.sin(rectangles[i].angle);
				double iV2x = -iV1y;
				double iV2y = iV1x;

				//change of basis
				double[][] newBasis = {
						{iV1x,iV2x},
						{iV1y,iV2y}
				};

				//the center of the circle is origo
				double[] corner1 = {rectangles[i].x1 - balls[j].x, rectangles[i].y1 - balls[j].y};
				double[] corner4 = {rectangles[i].x4 -  balls[j].x, rectangles[i].y4 - balls[j].y};

				//transformes the 2 diagonal corners of the rectangle to the new basis
				double[] transformedcorner1 = changeBasis(newBasis,corner1);
				double[] transformedcorner4 = changeBasis(newBasis,corner4);

				//calculates the closest point on the rectangle to a given ball
				//double xn = Math.max(transformedcorner1[0],Math.min(0,transformedcorner4[0]));
				//double yn = Math.max(transformedcorner1[1],Math.min(0,transformedcorner4[1]));

				double xn = clamp(0, transformedcorner1[0], transformedcorner4[0]);
				double yn = clamp(0, transformedcorner1[1], transformedcorner4[1]);

				int touple[] = new int[2];
				touple[0] = i;
				touple[1] = j;

				Coordinate collisionPoint = new Coordinate(xn, yn);
				collisionPoints.put(Arrays.toString(touple),collisionPoint);
			}
		}
		return collisionPoints;
	}

	private double clamp(double val, double min, double max) {
		return Math.max(min, Math.min(max, val));
	}

	private double[] changeBasis(double[][] basisMatrix, double[] vector) {
		double[][] inverseMatrix = invertMatrix(basisMatrix);
		return multiplyMatrixVector(inverseMatrix, vector);
	}

	private double[][] invertMatrix(double[][] basisMatrix) {
		double a = basisMatrix[0][0];
		double b = basisMatrix[0][1];
		double c = basisMatrix[1][0];
		double d = basisMatrix[1][1];

		double determinant = a * d - b * c;

		if (determinant == 0) {
			throw new IllegalArgumentException("Matrix is not invertible.");
		}

		double[][] inverse = new double[2][2];
		inverse[0][0] = d / determinant;
		inverse[0][1] = -b / determinant;
		inverse[1][0] = -c / determinant;
		inverse[1][1] = a / determinant;

		return inverse;
	}

	private double[] multiplyMatrixVector(double[][] inverseMatrix, double[] vector) {
		double[] result = new double[2];
		result[0] = inverseMatrix[0][0] * vector[0] + inverseMatrix[0][1] * vector[1];
		result[1] = inverseMatrix[1][0] * vector[0] + inverseMatrix[1][1] * vector[1];
		return result;
	}

	private void checkBalls(double deltaT) {


		HashMap< String, Double> ballPairMap = getBallDistances();

		for (int i = 0;  i<balls.length; i++) {
			for (int j = 0; j < balls.length; j++) {
				if(i != j) {
					int touple[] = new int[2];
					touple[0] = i;
					touple[1] = j;

					Ball oldBall1 = new Ball(balls[i]);
					Ball oldBall2 = new Ball(balls[j]);

					if (ballPairMap.get(Arrays.toString(touple)) < balls[i].radius + balls[j].radius) {

						// CollisionCenterLine is the line between the centerpoints of the balls.
						double collisionCenterLineX = balls[i].x - balls[j].x;
						double collisionCenterLineY = balls[i].y - balls[j].y;
						double collisionCenterLineVectorLength = Math.sqrt(Math.pow(collisionCenterLineX, 2) + Math.pow(collisionCenterLineY, 2));

						// Normalizes the vector
						collisionCenterLineX /= collisionCenterLineVectorLength;
						collisionCenterLineY /= collisionCenterLineVectorLength;

						// u1, u2, v1 and v2 refers to the velocity along the line between the balls center points
						// u1 is the speed of ball1 crashes into ball 2, reversed for u2 for ball2
						// v1 is the speed of ball1 after the collision, same applies to v2 for ball2
						// m1 is the mass for ball1, the same applies for ball2
						// mom is the momentum
						// rv is the relative velocity before the collision

						// Computes the dot product (procetion of the speed onto the line between the centerpoints).
						double u1 = (collisionCenterLineX * balls[i].vx) + (collisionCenterLineY * balls[i].vy);
						double u2 = (collisionCenterLineX * balls[j].vx) + (collisionCenterLineY * balls[j].vy);

						double m1 = balls[i].mass;
						double m2 = balls[j].mass;

						double mom = m1 * u1 + m2 * u2;
						double rv = u2 - u1;

						// these formulas for the new velocities after the collision while conserving the momentum and the kinetic energy.
						double v1 = (mom + m2 * rv) / (m1 + m2);
						double v2 = v1 - rv;

						// converts the velocities from being relative to the line between the balls center points to cooridinate system.
						balls[i].vx += collisionCenterLineX * (v1 - u1);
						balls[i].vy += collisionCenterLineY * (v1 - u1);
						balls[j].vx += collisionCenterLineX * (v2 - u2);
						balls[j].vy += collisionCenterLineY * (v2 - u2);

						// This method call will provide us with the distance between the balls if they continue to move with the same speed and direction for one more frame.
						double newDistance = nextBallDistance(deltaT,i,j);

						// This if-statement will hinder the balls from colliding a second time when they are moving away from each other after the collision but are still
						// overlapping. If the balls collide a second time, which sould not happend, then the balls will bounce towards each other.
						if (newDistance < ballPairMap.get(Arrays.toString(touple))) {
							// Switching back to the old speed and direction makes sure they will continue to move away from each other after the collision.
							balls[i] = oldBall1;
							balls[j] = oldBall2;
						}

					}
				}
			}
		}
	}

	private HashMap<String, Double> getBallDistances() {
		HashMap< String, Double> ballPairMap = new HashMap<>();
		for (int i = 0;  i<balls.length; i++){
			for (int j= 0; j< balls.length; j++){
				if(i != j) {

					double ballDistance = Math.sqrt(Math.pow(balls[i].x - balls[j].x, 2) + Math.pow(balls[i].y - balls[j].y, 2));
					int touple[] = new int[2];
					touple[0] = i;
					touple[1] = j;
					ballPairMap.put(Arrays.toString(touple),ballDistance);
				}
			}
		}
		return ballPairMap;
	}

	// Simulates the position of the balls if the continue to travel with the current speed and returns the distance between the balls.
	private double nextBallDistance(double deltaT,int i, int j) {
		double ball1nextX = balls[i].x + balls[i].vx * deltaT;
		double ball1nextY = balls[i].y + balls[i].vy * deltaT;
		double ball2nextX = balls[j].x + balls[j].vx * deltaT;
		double ball2nextY = balls[j].y + balls[j].vy * deltaT;
		return Math.sqrt(Math.pow(ball1nextX - ball2nextX, 2) + Math.pow(ball1nextY - ball2nextY, 2));
	}

	private double nextRectBallDistance(double deltaT, double centerX, double centerY, int i, int j) {
		double rectNextX = centerX + rectangles[i].vx * deltaT;
		double rectNextY = centerY + rectangles[i].vy * deltaT;
		double ballnextX = balls[j].x + balls[j].vx * deltaT;
		double ballnextY = balls[j].y + balls[j].vy * deltaT;
		return Math.sqrt(Math.pow((rectNextX - ballnextX),2) + Math.pow((rectNextY - ballnextY),2));
	}


	class Ball {
		
		Ball(double x, double y, double vx, double vy, double r, double m) {
			this.x = x;
			this.y = y;
			this.vx = vx;
			this.vy = vy;
			this.radius = r;
			this.mass = m;
		}

		// Copy constructor
		Ball(Ball ball) {
			this.x = ball.x;
			this.y = ball.y;
			this.vx = ball.vx;
			this.vy = ball.vy;
			this.mass = ball.mass;
			this.radius = ball.radius;
		}
		double x, y, vx, vy, radius, mass;
	}

	class Rectangle{

		//FÖRENKLA DETTA!!!!!!!!!!!!!!
		Rectangle(double px, double py, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4,
				  double vx, double vy, double rotation, double m, double a, double h, double w){
			this.positionX = px;
			this.positionY = py;
			this.x1 = x1;
			this.x2 = x2;
			this.y1 = y1;
			this.y2 = y2;
			this.x3 = x3;
			this.x4 = x4;
			this.y3 = y3;
			this.y4 = y4;
			this.vx = vx;
			this.vy = vy;
			this.rotation = rotation;
			this.mass = m;
			this.angle = a;
			this.height = h;
			this.width = w;

		}
		double positionX, positionY, x1, x2, x3, x4, y1, y2, y3, y4, vx, vy, rotation, mass, angle, height, width;

		//copy constructor
		Rectangle(Rectangle rectangle){
			this.positionX = rectangle.positionX;
			this.positionY = rectangle.positionY;
			this.x1 = rectangle.x1;
			this.x2 = rectangle.x2;
			this.y1 = rectangle.y1;
			this.y2 = rectangle.y2;
			this.x3 = rectangle.x3;
			this.x4 = rectangle.x4;
			this.y3 = rectangle.y3;
			this.y4 = rectangle.y4;
			this.vx = rectangle.vx;
			this.vy = rectangle.vy;
			this.rotation = rectangle.rotation;
			this.mass = rectangle.mass;
			this.angle = rectangle.angle;
			this.height = rectangle.height;
			this.width = rectangle.width;

		}
	}
}
