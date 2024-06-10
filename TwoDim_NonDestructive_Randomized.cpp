#include <stdio.h>
#include <ctime>
#include <stdlib.h>
#include <cmath>
#include <float.h>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <unistd.h>




long double Pi_Val = 3.14159265358979323846264338327950288419716939937510L;



std::string RandomString(const int len) {
    static const char alphanumerics[62] = {
            0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x41, 0x42,
            0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49, 0x4A, 0x4B, 0x4C, 0x4D, 0x4E,
            0x4F, 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5A,
            0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6A, 0x6B, 0x6C,
            0x6D, 0x6E, 0x6F, 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78,
            0x79, 0x7A
    };
    srand(time(0));

    std::string tmp;
    tmp.reserve(len);
    for (int i = 0; i < len; ++i) {
        tmp += alphanumerics[rand() % (sizeof(alphanumerics) - 1)];
    }
    return tmp;
}




typedef struct{
	long double x_axis;
	long double y_axis;
}Position;

typedef struct{
	long double x_axis;
	long double y_axis;
}TwoDVector;

typedef struct{
	long double x_displ;
	long double y_displ;
}CartesianDisplacement;

typedef struct{
	long double dist;
	long double angle;
}PolarDisplacement;

long double dotProduct(const Position& a, const Position& b){
    return a.x_axis * b.x_axis + a.y_axis * b.y_axis;
}

long double calculateDistance(const Position& a, const Position& b){
    return std::sqrt(((b.x_axis - a.x_axis) * (b.x_axis - a.x_axis)) + ((b.y_axis - a.y_axis) * (b.y_axis - a.y_axis)));
}

long double magnitudeSquared(const Position& p){
    return (p.x_axis * p.x_axis) + (p.y_axis * p.y_axis);
}


// Function to calculate distance from a point to a line segment
long double pointToSegmentDistance(const Position& p, const Position& start, const Position& end) {
    Position lineVec = {end.x_axis - start.x_axis, end.y_axis - start.y_axis};
    Position pointVec = {p.x_axis - start.x_axis, p.y_axis - start.y_axis};
    
    long double lineLenSquared = magnitudeSquared(lineVec);
    if (lineLenSquared == 0) return calculateDistance(p, start); // Line segment is a point

    // Consider the line extending the segment, parameterized as start + t (end - start).
    // We find projection of point p onto the line. 
    // It falls where t = [(p-start) . (end-start)] / |end-start|^2
    long double t = std::max(0.0L, std::min(1.0L, dotProduct(pointVec, lineVec) / lineLenSquared));
    Position projection = {start.x_axis + (t * lineVec.x_axis), start.y_axis + (t * lineVec.y_axis)}; // Projection falls on the segment
    return calculateDistance(p, projection);
}




CartesianDisplacement PolarToCartesian(PolarDisplacement Polardispl){
	CartesianDisplacement Neu;
	Neu.x_displ = 0.0L, Neu.y_displ = 0.0L;
	Neu.x_displ = Polardispl.dist * cosl(Polardispl.angle);
	Neu.y_displ = Polardispl.dist * sinl(Polardispl.angle);
	
	
	return Neu;
}

void ApplyCartFlight(Position* Position, CartesianDisplacement Flight){
	Position->x_axis = Position->x_axis + Flight.x_displ;
	Position->y_axis = Position->y_axis + Flight.y_displ;
	return;
}

void ApplyPolarFlight(Position* Position, PolarDisplacement PolarFlight){
	CartesianDisplacement CartesianFlight = PolarToCartesian(PolarFlight);
		
	Position->x_axis = Position->x_axis + CartesianFlight.x_displ;
	Position->y_axis = CartesianFlight.y_displ + Position->y_axis;
	return;
}







long double PowerDistribution(long double tailpower, long double SightLimit){
	

	std::mt19937_64 generator(rand());
	std::uniform_real_distribution<long double> distribution(0.0L,1.0L); //Standard Uniform
	
	long double UniformRandNumb = distribution(generator);
	long double x = (SightLimit * powl(1.0L - UniformRandNumb, -1.0L / (tailpower - 1.0L)));
	
	
	//std::cout << "Polar Dist: " << x << std::endl;
	return x;
}


long double TurnAngle(){
	
	std::mt19937_64 generator(rand());
	std::uniform_real_distribution<long double> distribution(0.0L, 131072.0L * Pi_Val);
	
	long double alfa = distribution(generator);
	//std::cout << "Alfa: " << alfa << std::endl;
	return alfa;
}



PolarDisplacement Flight(long double tailpower, long double SightLimit){
	if(tailpower > 0.0L){
		if(((!(tailpower > 1.000000L)) || (tailpower > 3.000000L))){
			printf("Wrong tail value");	
		}
	}
	
	
	PolarDisplacement NeuFlight;
	NeuFlight.dist = PowerDistribution(tailpower, SightLimit);
	NeuFlight.angle = TurnAngle();
	
	return NeuFlight;
}




std::vector<Position> findCircleLineIntersections(Position circleCenter, long double radius, Position lineStart, Position lineEnd){
	std::vector<Position> intersections;
	//std::cout << "Node Position: " << circleCenter.x_axis << " " << circleCenter.y_axis << std::endl;
	

    long double dx = lineEnd.x_axis - lineStart.x_axis;
    long double dy = lineEnd.y_axis - lineStart.y_axis;
    long double A = (dx * dx) + (dy * dy);
    long double B = 2 * (dx * (lineStart.x_axis - circleCenter.x_axis) + dy * (lineStart.y_axis - circleCenter.y_axis));
    long double C = (lineStart.x_axis - circleCenter.x_axis) * (lineStart.x_axis - circleCenter.x_axis) +
                    (lineStart.y_axis - circleCenter.y_axis) * (lineStart.y_axis - circleCenter.y_axis) -
                    radius * radius;

    long double det = B * B - 4 * A * C;
    if (det < 0.00000000000000000000000000000000L) {
        // No real solutions.
        return intersections;
    } else if (det < 0.000000000001L) { //Chosen Tolerance limit of 0
        // One solution.
        //std::cout << "One Solution for: " << circleCenter.x_axis << " " << circleCenter.y_axis << std::endl;
        //std::cout << "Start: " << lineStart.x_axis << " " << lineStart.y_axis << std::endl;
        //std::cout << "End: " << lineEnd.x_axis << " " << lineEnd.y_axis << std::endl;
        long double t = -B / (2 * A);
        intersections.push_back({lineStart.x_axis + (t * dx), lineStart.y_axis + (t * dy)});
    } else {
        // Two solutions.
        //std::cout << "Double Solution for: " << circleCenter.x_axis << " " << circleCenter.y_axis << std::endl;
        //std::cout << "Start: " << lineStart.x_axis << " " << lineStart.y_axis << std::endl;
        //std::cout << "End: " << lineEnd.x_axis << " " << lineEnd.y_axis << std::endl;
        long double t = (-B + std::sqrt(det)) / (2 * A);
        long double t2 = (-B - std::sqrt(det)) / (2 * A);
        Position intersection1 = {lineStart.x_axis + (t * dx), lineStart.y_axis + (t * dy)};
        Position intersection2 = {lineStart.x_axis + (t2 * dx), lineStart.y_axis + (t2 * dy)};
        if (fabs(calculateDistance(lineStart, intersection1) + calculateDistance(intersection1, lineEnd) - calculateDistance(lineStart, lineEnd)) < 0.001L) {
            intersections.push_back(intersection1);
        }
        if (fabs(calculateDistance(lineStart, intersection2) + calculateDistance(intersection2, lineEnd) - calculateDistance(lineStart, lineEnd)) < 0.001L) {
            intersections.push_back(intersection2);
        }
    }

    return intersections;
}





Position findEncounteredNode(const std::vector<Position>& nodes, const Position start, PolarDisplacement PD, long double VisionRange, size_t* foundflag){

    if(foundflag == NULL){
		
    	// Calculate end point of the forager's movement
    	Position end = start;
    	ApplyPolarFlight(&end, PD);
    	//std::cout << "Future Position: " << end.x_axis << " " << end.y_axis << std::endl;
    
    
    	Position closestEncounter;
    	long double minDistance = PD.dist + 0.000001;
    	//std::cout << "Size: " << nodes.size() << std::endl;

    	for (const auto& node : nodes){
        	auto intersections = findCircleLineIntersections(node, VisionRange, start, end);

        	for (const auto& point : intersections){
        		std::cout << "intersections: " << point.x_axis << " " << point.y_axis << std::endl;
            	long double distance = calculateDistance(start, point);
            	if(distance < minDistance) {
                	minDistance = distance;
                	closestEncounter = point;
            	}
        	}
    	}

    	return closestEncounter; //ignored if no encounter
	}
	
	*foundflag = -1; //Ensure flag is reset
	
    // Calculate end point of the forager's movement
    Position end = start;
    ApplyPolarFlight(&end, PD);
    std::cout << "Future Position: " << end.x_axis << " " << end.y_axis << std::endl;
    
    
    Position closestEncounter;
    long double minDistance = PD.dist + 0.00000000001; //Ensure intersection points are within displacement range


    for (size_t i = 0; i < nodes.size(); i++){
    	//std::cout << "Tested: " << nodes[i].x_axis << " " << nodes[i].y_axis << std::endl;
        auto intersections = findCircleLineIntersections(nodes[i], VisionRange, start, end);

        for (const auto& point : intersections){
        	//std::cout << "intersections: " << point.x_axis << " " << point.y_axis << std::endl;
            long double distance = calculateDistance(start, point);
            if (distance < minDistance){  //Update better fit candidate
            	//std::cout << "Fitter intersection: " << point.x_axis << " " << point.y_axis << std::endl;
            	*foundflag = i;  //Save index of success
                minDistance = distance;
                closestEncounter = point;
            }
        }
    }

    return closestEncounter; //ignored if no encounter
}









int main(void){
	srand (time(NULL));
	std::mt19937_64 generatorDos(rand());
	long double MapBounds = 10000.0L;
    std::uniform_real_distribution<long double> Unifff(-1.0L * MapBounds, MapBounds);
    std::uniform_real_distribution<long double> UnifffInit(-1500.0L, 1500.0L);
	std::uniform_real_distribution<long double> UnifffTail(1.02L, 2.98L);             //Pick \mu range
	
	
	Position PositionForager;
	PositionForager.x_axis = 0.0L;
	PositionForager.y_axis = 0.0L;
	
	
	
	
	long double PowerTail = UnifffTail(generatorDos);       //randomly chosen
	long double MaxSight =  0.25L;
	
	std::ofstream fout;
	std::string FileName;
	if(PowerTail > 3.0L){
		FileName = "NonDestructiveBrownianMot_" + RandomString(16);
	}else{
		FileName = "NonDestructiveLevyFlights_" + RandomString(16);
	}
	
	fout.open(FileName.c_str());
	
	std::cout << "Non Destructive Forager " << FileName.c_str() << std::endl;
	std::cout << PowerTail << " tail" << std::endl;
	std::cout << MaxSight << " Sight Radius" << std::endl;
	fout << PowerTail << " tail" << std::endl;
	fout << MaxSight << " Sight Radius" << std::endl;

	
	//Initialize TargetGrid
	//Let's consider target area from -500 to +500 in both axis -> 1000 * 1000 = 1000000 area units
	const unsigned int nrolls = 8000000; // total number of nodes
	std::vector<Position> NodePositions; 
	NodePositions.reserve(nrolls + 1);  
	
	

    
    std::cout << "Initializing Nodes with density " << nrolls/(2.0L * MapBounds * 2.0L * MapBounds) << std::endl;
	
	std::cout << "Mean Free Path " << 1/(2 * MaxSight * ((nrolls/(2.0L * MapBounds * 2.0L * MapBounds)))) << std::endl;
	
		
	for(unsigned int i = 0; i < nrolls; i++){
		Position NodePositionA;
		NodePositionA.x_axis = Unifff(generatorDos);
		NodePositionA.y_axis = Unifff(generatorDos);
		//NodePositions[i] = NodePositionA;
		NodePositions.push_back(NodePositionA);
		//std::cout << "Node Pos: " << NodePositions[i].x_axis << " " << NodePositions[i].y_axis << std::endl;
		
	}
	
	std::cout << "Nodes Initialized" << std::endl;
	
	
	size_t Foundflag = -1;
	
	unsigned int MaxIter = 1000;
	unsigned int totalflights = 0, itter = 0;
	long double totalDispl = 0.0L;
		do{
		
		std::cout << std::endl;
		do{ //Scoped
			unsigned int localflights = 0;
			long double LocalDispl = 0.0L;
			std::vector<Position> NonDestructivePlaceholderNodePositions; 
			NonDestructivePlaceholderNodePositions.reserve(nrolls + 1);  
			for (size_t i=0; i<NodePositions.size(); i++){    ////TODO: This could be simplified by using plain NodePositions
			    NonDestructivePlaceholderNodePositions.push_back(NodePositions[i]); 
			}
			Position RecentlyCollectedNode;
			int RecentlyCollected = -1;
			do{
				std::cout << "Start Pos: " << PositionForager.x_axis << " " << PositionForager.y_axis << std::endl;
				if((PositionForager.x_axis < MapBounds) && (PositionForager.x_axis > (-1.0L * MapBounds)) && (PositionForager.y_axis < MapBounds) && (PositionForager.y_axis > (-1.0L * MapBounds))){
					//Forager still bounded
				}else{
					std::cout << "Forager outside bounds, halting..." << std::endl;
					fout << std::endl;
					fout << "Forager outside bounds, halting..." << std::endl;
					fout << std::endl;
					exit(42);
					break;
				}
				localflights++;
				PolarDisplacement AFlight = Flight(PowerTail, MaxSight);
				std::cout << "Dist: " << AFlight.dist << std::endl;
				auto FirstEncounter = findEncounteredNode(NonDestructivePlaceholderNodePositions, PositionForager, AFlight, MaxSight, &Foundflag);
				if(Foundflag != (size_t)-1){
					if(RecentlyCollected != -1){  //Restore previous collected node to the back position, so that no index errors occur
						RecentlyCollected = -1;
						NonDestructivePlaceholderNodePositions.push_back(RecentlyCollectedNode);
					}
					std::cout << "Encounter: " << FirstEncounter.x_axis << " " << FirstEncounter.y_axis << std::endl;  //This contains the point at which travel line intersects with MaxSight influence circle
					std::cout << "Node: " << NonDestructivePlaceholderNodePositions[Foundflag].x_axis << " " << NonDestructivePlaceholderNodePositions[Foundflag].y_axis << std::endl;  //This contains the point at which travel line intersects with MaxSight influence circle
					auto DirectDist = calculateDistance(NonDestructivePlaceholderNodePositions[Foundflag], PositionForager);
					if(DirectDist < MaxSight){ 
						//Node already within sight
						LocalDispl += DirectDist;
					}else{
						//Node outside sight
						LocalDispl += (MaxSight + calculateDistance(FirstEncounter, PositionForager));
					}
					PositionForager = NonDestructivePlaceholderNodePositions[Foundflag];  //Update position
					itter++;   //Success increase
					RecentlyCollected = 1;  //Update Collection Flag
					RecentlyCollectedNode = NonDestructivePlaceholderNodePositions[Foundflag];  //Save erased node for future restoring
					NonDestructivePlaceholderNodePositions.erase(NonDestructivePlaceholderNodePositions.begin() + Foundflag);  //Remove foraged node from next immediate search
					Foundflag = -1;	//Update encounter flag
					
					if((localflights%100) == 50){  //Scoped, save results every some successes
						fout << std::endl;
						fout << std::endl;
						fout << "Success Backup: " << std::endl;
						fout << localflights << " flights" << std::endl;
						fout << itter << " successes" << std::endl;
						fout << "Traversed " << LocalDispl << " Distance Units" << std::endl;
  	
  						double AverageFlightCount = (1.0 * localflights * 1.0)/(1.0 * itter * 1.0);
  						long double AverageSuccessLength = LocalDispl/(1.0L * itter * 1.0L);
  						long double AverageFlightLength = LocalDispl/(1.0L * localflights * 1.0L);
  	
  						fout << "Current Averages: " << std::endl;
						fout << AverageFlightCount << " flights per success" << std::endl;
						fout << AverageSuccessLength << " Distance Units per success" << std::endl;
  						fout << AverageFlightLength << " Distance Units per flight" << std::endl;
						fout << std::endl;
						fout << std::endl;
					}
									
				}else{
					if(RecentlyCollected != -1){  //Restore node previously collected and update flag 
						RecentlyCollected = -1;
						NonDestructivePlaceholderNodePositions.push_back(RecentlyCollectedNode);
					}
					ApplyPolarFlight(&PositionForager, AFlight);
					LocalDispl += AFlight.dist;
				}
				
				/*
				//Keep track of progress
				if((localflights%50) == 49){
					fout << std::endl;
					fout << std::endl;
					fout << "Treshold Backup: " << localflights << " " << std::endl;
					fout << itter << " successes" << std::endl;
					fout << "Traversed " << LocalDispl << " Distance Units" << std::endl;
  	
  					double AverageFlightCount = (1.0 * localflights * 1.0)/(1.0 * itter * 1.0);
  					long double AverageSuccessLength = LocalDispl/(1.0L * itter * 1.0L);
  					long double AverageFlightLength = LocalDispl/(1.0L * localflights * 1.0L);
  	
  					fout << "Current Averages: " << std::endl;
					fout << AverageFlightCount << " flights per success" << std::endl;
					fout << AverageSuccessLength << " Distance Units per success" << std::endl;
  					fout << AverageFlightLength << " Distance Units per flight" << std::endl;
					fout << std::endl;
					fout << std::endl;
				}
				*/
				
				//std::cout << "New Pos: " << PositionForager.x_axis << " " << PositionForager.y_axis << std::endl;
			}while((localflights < 1000 * MaxIter) && (itter < MaxIter)); //Halt at MaxIter successes or at the total of 1000 * MaxIter decisions
			
			if(localflights > 0){
				totalDispl += LocalDispl;
				totalflights += localflights;
			}
			
		}while(false);
		
		if(!(totalDispl < (1000 * MaxIter))){
			std::cout << "Forager stuck for too long, halting..." << std::endl;
			fout << std::endl;
			fout << "Forager stuck for too long, halting..." << std::endl;
			fout << std::endl;
			break;
		}
		
	}while(itter < MaxIter);

	std::cout  << std::endl;
	std::cout << "For: " << itter << " successes" << std::endl;
	std::cout << "it took: " << totalflights << " flights" << std::endl;
	std::cout << "For a total of " << totalDispl << " Distance Units Traversed" << std::endl;
	std::cout  << std::endl;	
  	
  	double AverageFlightCount = (1.0 * totalflights * 1.0)/(1.0 * itter * 1.0);
  	long double AverageSuccessLength = totalDispl/(1.0L * itter * 1.0L);
  	long double AverageFlightLength = totalDispl/(1.0L * totalflights * 1.0L);
  	
  	std::cout << "For an average of: " << std::endl;
	std::cout << AverageFlightCount << " flights" << std::endl;
	std::cout << AverageSuccessLength << " Distance Units per success" << std::endl;
  	std::cout << AverageFlightLength << " Distance Units per flight" << std::endl;
	
	fout << std::endl;
	fout << std::endl;
	fout << std::endl;
	fout << std::endl;
	
	
	fout  << std::endl;
	fout << PowerTail << " tail" << std::endl;
	fout << MaxSight << " Sight Radius" << std::endl;
	fout << "For: " << itter << " successes" << std::endl;
	fout << "it took: " << totalflights << " flights" << std::endl;
	fout << "For a total of " << totalDispl << " Distance Units Traversed" << std::endl;
	fout << "For an average of: " << std::endl;
	fout << AverageFlightCount << " flights" << std::endl;
	fout << AverageSuccessLength << " Distance Units per success" << std::endl;
  	fout << AverageFlightLength << " Distance Units per flight" << std::endl;
  	
  	long double EtaEfficiency = 1.0L/(AverageFlightCount * AverageFlightLength);
  	fout << EtaEfficiency << " Eta Efficiency" << std::endl;




	fout  << std::endl;	
	fout.close();
	
	return 0;
}
