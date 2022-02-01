import java.util.*;
import java.lang.Math;

class City{
	public double x;
	public double y;
	
	
	public void setx(double x) {
		this.x = x;
	}
	
	public void sety(double y) {
		this.y = y;
	}
	
	@Override
	public String toString() {
		return "(" + x + "," + y + ")";
	}

}


public class Tour {
	
	public City [] cities;

	Tour (City[] cities){
		this.cities = cities;
	}
	
	public static double distance (City c1, City c2) {
		return Math.sqrt(Math.pow((c2.x - c1.x),2) + Math.pow((c2.y - c1.y),2));
	}
	
	
	public static double cost (City [] cities) {
		double total = 0;
		for (int i = 0; i<cities.length-1; i++) {
			total += distance(cities[i], cities[i+1]);
		}
		total += distance(cities[cities.length - 1], cities[0]); //loop back
		return total;
	}
	
	//TSP represented as an instance of Tour
	public City[] TSP (City [] cities, Random rand) {
		int i = 0;
		for (City city:cities) {
			double r1 = rand.nextDouble();
			double r2 = rand.nextDouble();
			city.setx(r1);
			city.sety(r2);
			cities[i] = city;
			i++;
		}
		
		return cities;
	}
	
	//Generate all TSP instances 
	public City[][] generateAllTSP(int size, City [] cities) {
		
		City [][] c = new City [size][];
		for (int i = 0; i < size; i++) {
			Random rand = new Random();
			c[i] = TSP(cities, rand);	
		}
		return c;
	}
	
	
	public static double [][] computeDistance (City [] cities){
		
		double [][] result = new double [cities.length][cities.length];
		for (int i = 0; i < cities.length; i++) {
			for (int j = 0; j < cities.length; j++) {
				result[i][j] = distance(cities[i], cities[j]);
			}
		}
		return result;
	}
	
	public double bruteForce (City[] cities, double [][] distances, boolean[] visited, double cycleCost, double cost, int position, int iteration) {
		//Memoization
		if (iteration == cities.length && distances[position][0] > 0) {
			cycleCost = Math.min(cycleCost, cost + distances[position][0]);
			return cycleCost;
		}
		
		for (int i = 0; i<distances.length; i++) {
			if (visited[i] == false && distances[position][i] > 0) {//make sure of zero distance (same city)
				visited[i] = true; //visit city
				cycleCost = bruteForce(cities, distances, visited, cycleCost, cost + distances[position][i], i, (iteration + 1));
				visited[i] = false; //unvisit city
			}
		}
		return cycleCost;
	}
	
	public List<List<City>> generateNeighbors(City[] c){
		
		List<City> cList = new ArrayList<City>();
		cList = Arrays.asList(c);
		List<List<City>> neighbors = new ArrayList<List<City>>();

		for (int i = 0; i<c.length; i++) {
			for(int j = 0; j<c.length;j++) {
				if (j == i) {continue;}
				
				List<City> copy = new ArrayList<City>();
				copy.addAll(cList); 
				City ci = copy.get(i);
				copy.set(i, copy.get(j));
				copy.set(j, ci);
				
				if (!neighbors.contains(copy)) {
					neighbors.add(copy);
				}
			}
		}
		return neighbors;
	}
	
	public City [] change (List<City> c) {
		City [] cities = new City [c.size()]; 
		for (int i = 0; i < c.size(); i++) {
			cities[i] = c.get(i);
		}
		return cities;
	}
	
	public List<City> getBestNeighbor(List<List<City>> neighbors){
		double best = Double.MAX_VALUE;
		List<City> bestCities = new ArrayList<City>();
		for (int i = 0; i < neighbors.size(); i++) {//iterate through all neighbors
			City [] cities1 = change(neighbors.get(i));
			if (cost(cities1) < best) {
				best = cost(cities1); 
				bestCities = Arrays.asList(cities1);
			}
		}
		return bestCities;
	}
	
	
	public double hillClimbing(List<City> c) {
		City [] cities = change(c);
		double distance = cost(cities);
		List<List<City>> neighbors = generateNeighbors(cities);
		List<City> bestCities = getBestNeighbor(neighbors);
		double best = cost(change(bestCities));
		
		while (best < distance) {
			distance = best;
			neighbors = generateNeighbors(change(bestCities));
			bestCities = getBestNeighbor(neighbors);
			best = cost(change(bestCities));
		}
		return distance;
	}
	
	public void solveTSP(City [] c) {
		//part a variables
		double min = Double.MAX_VALUE;
		double max = -1.0;
		double sum = 0.0;
		double sd = 0.0;
		double [] tourDistances = new double [100];
		
		//part b variables
		double randomMin = Double.MAX_VALUE;
		double randomMax = -1.0;
		double randomSum = 0.0;
		double randomSD = 0.0;
		double [] randomTourDistances = new double [100];
		int counter = 0;
		
		//part c
		double hillMin = Double.MAX_VALUE;
		double hillMax = -1.0;
		double hillSum = 0.0;
		double hillSD = 0.0;
		double [] hillDistances = new double [100];
		int hillCounter = 0;
		
		
		for (int i = 0; i<100; i++) {//for each tour instance
			
			//part a
			City [] cities = generateAllTSP(100, c)[i];
			boolean [] visited = new boolean[cities.length];
			visited[0] = true;
			double [][] distances = computeDistance(cities);
			double distance = Double.MAX_VALUE;
			distance = bruteForce(cities, distances, visited, distance, 0.0, 0, 1);
			distance = (double) Math.round(distance * 1000000d) / 1000000d;
			
			if (distance < min) {
				min = distance;
			}
			if (distance > max) {
				max = distance;
			}
			sum += distance;
			tourDistances[i] = distance;
			
			//Part b
			List<City> citiesList = Arrays.asList(cities);
			Collections.shuffle(citiesList);
			citiesList.toArray(cities);
			double randomDist = cost(cities);
			randomDist = (double) Math.round(randomDist * 1000000d) / 1000000d;
			
			if (randomDist == distance) {
				counter++;
			}
			if (randomDist < randomMin) {
				randomMin = randomDist;
			}
			if (randomDist > randomMax) {
				randomMax = randomDist;
			}
			randomTourDistances[i] = randomDist;
			randomSum += randomDist;	
			
			//part c
			double hillDistance = hillClimbing(citiesList);
			hillDistance = (double) Math.round(hillDistance * 1000000d) / 1000000d;
			
			if (hillDistance < hillMin) {
				hillMin = hillDistance;
			}
			if (hillDistance > hillMax) {
				hillMax = hillDistance;
			}
			if (hillDistance == distance) {
				hillCounter++;
			}
			hillDistances[i] = hillDistance;
			hillSum += hillDistance;
		}
		
		double mean = sum/100;
		double randomMean = randomSum/100;
		double hillMean = hillSum/100;

		for (int i = 0; i<100; i++) {
			sd += Math.abs(tourDistances[i] - mean);
			randomSD += Math.abs(randomTourDistances[i] - randomMean);
			hillSD += Math.abs(hillDistances[i] - hillMean);
		}
		
		sd = sd/10;
		randomSD = randomSD/10;
		hillSD = hillSD/10;

		//print out result
		System.out.println("Part a: ");
		System.out.println("Mean: " + mean);
		System.out.println("Min: " + min);
		System.out.println("Max: " + max);
		System.out.println("Standard deviation: " + sd);
		//System.out.println(Arrays.toString(tourDistances));
		System.out.println();
		System.out.println("Part b: ");
		System.out.println("Mean: " + randomMean);
		System.out.println("Min: " + randomMin);
		System.out.println("Max: " + randomMax);
		System.out.println("Standard deviation: " + randomSD);
		System.out.println("Optimal solution found: " + counter);
		//System.out.println(Arrays.toString(randomTourDistances));
		System.out.println();
		System.out.println("Part c: ");
		System.out.println("Mean: " + hillMean);
		System.out.println("Min: " + hillMin);
		System.out.println("Max: " + hillMax);
		System.out.println("Standard deviation: " + hillSD);
		System.out.println("Optimal solution found: " + hillCounter);
		//System.out.println(Arrays.toString(hillDistances));
	}
	
	
	public void solveTSP1(City[] c) {
		//part b variables
		double randomMin = Double.MAX_VALUE;
		double randomMax = -1.0;
		double randomSum = 0.0;
		double randomSD = 0.0;
		double [] randomTourDistances = new double [100];
		
		//part c
		double hillMin = Double.MAX_VALUE;
		double hillMax = -1.0;
		double hillSum = 0.0;
		double hillSD = 0.0;
		double [] hillDistances = new double [100];
		
		for (int i = 0; i<100; i++) {//for each tour instance
			
			//Part b
			City [] cities = generateAllTSP(100, c)[i];
			List<City> citiesList = Arrays.asList(cities);
			Collections.shuffle(citiesList);
			citiesList.toArray(cities);
			double randomDist = cost(cities);
			randomDist = (double) Math.round(randomDist * 1000000d) / 1000000d;
			
			if (randomDist < randomMin) {
				randomMin = randomDist;
			}
			if (randomDist > randomMax) {
				randomMax = randomDist;
			}
			randomTourDistances[i] = randomDist;
			randomSum += randomDist;	
			
			//part c
			double hillDistance = hillClimbing(citiesList);
			hillDistance = (double) Math.round(hillDistance * 1000000d) / 1000000d;
			
			if (hillDistance < hillMin) {
				hillMin = hillDistance;
			}
			if (hillDistance > hillMax) {
				hillMax = hillDistance;
			}
			hillDistances[i] = hillDistance;
			hillSum += hillDistance;
		}
		double randomMean = randomSum/100;
		double hillMean = hillSum/100;

		for (int i = 0; i<100; i++) {
			randomSD += Math.abs(randomTourDistances[i] - randomMean);
			hillSD += Math.abs(hillDistances[i] - hillMean);
		}
		
		randomSD = randomSD/10;
		hillSD = hillSD/10;

		//print out result
		System.out.println();
		System.out.println("100 cities: ");
		System.out.println("Part b: ");
		System.out.println("Mean: " + randomMean);
		System.out.println("Min: " + randomMin);
		System.out.println("Max: " + randomMax);
		System.out.println("Standard deviation: " + randomSD);
		//System.out.println(Arrays.toString(randomTourDistances));
		System.out.println();
		System.out.println("Part c: ");
		System.out.println("Mean: " + hillMean);
		System.out.println("Min: " + hillMin);
		System.out.println("Max: " + hillMax);
		System.out.println("Standard deviation: " + hillSD);
		//System.out.println(Arrays.toString(hillDistances));
	}
	
	public static void main(String[] args) {
		//generate 7 cities
		City [] c = new City[7];
		for (int i = 0; i < 7; i++) {
			City c1 = new City();
			c[i] = c1;
		}
		//generate 100 cities
		City [] c100 = new City[100];
		for (int j = 0; j < 100; j++) {
			City c1 = new City();
			c100[j] = c1;
		}
		
		//create tour
		Tour t = new Tour(c);
		t.solveTSP(c);
		Tour t100 = new Tour(c100);
		t100.solveTSP1(c100);
	}

}
