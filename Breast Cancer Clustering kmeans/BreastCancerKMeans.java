package BreastCancerKMeans;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
public class BreastCancerKMeans {
	
	Data[] cancerData;
	
	public  int k;
	public int l;
	public Double iterNum;
	Random random = new Random();
		
	ArrayList<Cluster[]> centroidArray ;
	
	public BreastCancerKMeans(String file,int k1, int l1){
		BufferedReader reader;
		String line;
		String[] splitLine;
		Double[] attrValues;
		int recordCount;
		k=k1;
		l=l1;
		iterNum=0.00001;
		centroidArray = new ArrayList<Cluster[]>();
		try {
			
			reader = new BufferedReader(new FileReader(file));
		      
		     // Count the number of records in csv file
		    for(recordCount = 0; reader.readLine( ) != null; recordCount++);
		      
		      // Close and then re-open the file now that we know its length
		    reader.close();
		      
			reader = new BufferedReader(new FileReader(file));
		    
			cancerData = new Data[recordCount];
			attrValues = null;
		     
			for(int i=0; i< recordCount; i++){ 
				line=reader.readLine();
				splitLine = line.split( "," );
				
				// The first entry in the parts array is the scn, the rest
		        // are attribute values of the instance
		        if((attrValues == null) || (attrValues.length != 
		                                   (splitLine.length - 1)))
		        	attrValues = new Double[splitLine.length - 2];
		        for(int j = 0; j < attrValues.length; j++)
		        	attrValues[j] = Double.parseDouble(splitLine[j + 1]);
		        
		        // Finally, create the Gene in the array
		        cancerData[i] = new Data(splitLine[0], attrValues);
		        cancerData[i].setClassLabel(Integer.parseInt(splitLine[splitLine.length - 1]));
			}
			
			 // close the file
		    reader.close();
		      
		} catch (IOException e) {
			e.printStackTrace();
		}	
		
	}
	
	public Data[] getData() 
	  {
	    return cancerData; 
	  }
	
	public ArrayList<Cluster[]> getClusters() 
	  {
		
	    return centroidArray;
	  }
	  
	
	public void getInitialCentroids(){
		//Initial random centroid points from domain of data
		Cluster[] initialRandCentroids1 = new Cluster[k];
		initialRandCentroids1=getKRandomPoints();
		
		try {
			centroidArray.set(0, initialRandCentroids1) ;
		} catch ( IndexOutOfBoundsException e ) {
			centroidArray.add(0, initialRandCentroids1);
		}
		
		//Initial random centroid points from domain of data
		Cluster[] initialRandCentroids = new Cluster[k];
		initialRandCentroids=getKRandomPoints();
		
		//Assign data points to initial random centroids
		initialRandCentroids = assignDataToMultipleCentroids(initialRandCentroids);
		
		
		try {
			centroidArray.set(1, initialRandCentroids) ;
		} catch ( IndexOutOfBoundsException e ) {
			centroidArray.add(1, initialRandCentroids);
		}
		
	}
		
	
	public Cluster[] getKRandomPoints() {
		List<Data> kRandomPoints = new ArrayList<Data>();

		boolean[] alreadyChosen = new boolean[cancerData.length];
		int size = cancerData.length;
		Cluster[] clusters = new Cluster[k];
		for (int i = 0; i < k; i++) {
			int index = -1, r = random.nextInt(size--) + 1;
			for (int j = 0; j < r; j++) {
				index++;
				while (alreadyChosen[index])
					index++;
			}
			kRandomPoints.add(cancerData[index]);
			
			clusters[i] = new Cluster();
			clusters[i].assignCentroid(cancerData[index]);
			
		   alreadyChosen[index] = true;
		  }
				  return clusters;
		 }
		  
	
	
	public double calculateDistance(Double[] attrValues, Double[] doubles)
    {
        double Sum = 0.0;
        for(int i=0;i<attrValues.length;i++) {
           Sum = Sum + Math.pow((attrValues[i]-doubles[i]),2.0);
        }
        return Math.sqrt(Sum);
    }
	
	public Cluster[] assignDataCentroids(Cluster[] initialRandCentroids1){
		Double[] d;
		int min;
		for(int i=0;i<cancerData.length;i++){
			Double[] attrValues = cancerData[i].getAttributeValues();
			min = 0;
			d = new Double[initialRandCentroids1.length];
			for(int j=0;j<initialRandCentroids1.length;j++){
				
				d[j] = calculateDistance(attrValues,initialRandCentroids1[j].getCentroid().getAttributeValues());
				if(j!=0){
					if(d[j]<d[j-1]){
						min = j;
					}
				}
				
			}
			initialRandCentroids1[min].addElement(cancerData[i]);
		}
		return initialRandCentroids1;
		
	}
	
	public Cluster[] assignDataToMultipleCentroids(Cluster[] initialRandCentroids2){
		Double[] d;
		int min;
		
		
		for(int m=1; m<l+1; m++){
			for(int i=0;i<cancerData.length;i++){
				ArrayList<Integer> minDistanceArray = new ArrayList<>();
				Double[] attrValues = cancerData[i].getAttributeValues();
				min = 0;
				d = new Double[initialRandCentroids2.length];
				if(m==1){
					
				}else{
					minDistanceArray = (cancerData[i].getClusterIndex());
				}
				for(int j=0;j<initialRandCentroids2.length;j++){
					
					if(minDistanceArray.contains(j)){
						d[j]=999999999d;
					}else{
						
						d[j] = calculateDistance(attrValues,initialRandCentroids2[j].getCentroid().getAttributeValues());
						if(j!=0){
							if(d[j]<d[min]){
								min = j;
							}
						}
					}
					
					
				}
				cancerData[i].setclusterIndex(min, m-1);
				initialRandCentroids2[min].addElement(cancerData[i]);
				
			}
		}
		return initialRandCentroids2;
	}
	public void updateCentroid(){
		int centroidArrayLength; 
		Double interCentroidDist;
		Double flag = 0.0;
		Double sum;
		Double sumCurrent = 0.0;
		int iter = centroidArray.size();
		Double sumPrevCentroid = 0.0;
		do{
			
			Cluster[] clusterTemp = new Cluster[k];
			for(int i=0;i<k;i++){
				Cluster cluster1 = centroidArray.get(iter-1)[i];
				Data dataTemp = cluster1.getMean();
				clusterTemp[i] = new Cluster();
				clusterTemp[i].assignCentroid(dataTemp);
				
			}
			//Assign data points to updated Centroids
			clusterTemp = assignDataToMultipleCentroids(clusterTemp);
						
			try {
				centroidArray.set(iter, clusterTemp) ;
			} catch ( IndexOutOfBoundsException e ) {
				centroidArray.add(iter, clusterTemp);
			}
			
			centroidArrayLength = centroidArray.size();
			
			if (centroidArrayLength > 2){
				
				Double[] interCentroidDistance = new Double[2];
					for(int z=0;z<2;z++){
						sumCurrent = 0.0;
						Double [] minIntraCentroidDistance = new Double[k];
						for(int n=0;n<k;n++){
							Double[] intraCentroidDistance = new Double[k];
							for(int p=0;p<k;p++){
								sumCurrent = sumCurrent + calculateDistance((centroidArray.get(centroidArrayLength-1-z))[n].getCentroid().getAttributeValues(),(centroidArray.get(centroidArrayLength-2-z))[p].getCentroid().getAttributeValues());
							}
														
						}
						
						interCentroidDistance[z] = sumCurrent;
												
					}
					
					flag = Math.abs(interCentroidDistance[0] - interCentroidDistance[1]);
					sumPrevCentroid  = interCentroidDistance[1];
					System.out.println("Sum of current centroids:"+interCentroidDistance[0]);
					System.out.println("Sum of prev centroids:"+interCentroidDistance[1]);
					iter++;
			}
		}while(flag > ((iterNum)*(sumPrevCentroid)));
		
	}
	
	public static void main(String[] args) {
		int k=2;
		int l2=1;
		String filePath = "";
		
		BreastCancerKMeans breastCancerKmeans = new BreastCancerKMeans(filePath,k,l2);
		Double[] attrValues;
		Data[] cancerData = breastCancerKmeans.getData();
		breastCancerKmeans.getInitialCentroids();
		breastCancerKmeans.updateCentroid();
			
		ArrayList<Cluster[]> clustersData=breastCancerKmeans.getClusters();
		HashSet<Data> clusterValues = new HashSet<Data>();
		Double[] clusterValues1;
		Double[] centroidValues;
			
		
		int[] benignCount = new int[2];
		int[] malignCount = new int[2];
		
		double[] error = new double[2];
			Cluster[] l = new Cluster[k];
			l = clustersData.get(clustersData.size()-1);
			for(int i=0;i<l.length;i++){
				System.out.println("Centroid scn:"+i+""+l[i].getCentroid().getscn());
				System.out.println("***********************************************");
				System.out.println("Centroid "+i+" attribute values:");
				centroidValues = l[i].getCentroid().getAttributeValues();
				System.out.print("("+centroidValues[0]);
				for(int r=1;r<centroidValues.length;r++){
					System.out.print(","+centroidValues[r]);
				}
				System.out.println(")");
				System.out.println("***********************************************");
				System.out.println("Centroid cluster size:"+l[i].getClusterElements().size());
				
				for(Data d:l[i].getClusterElements()){
				
					clusterValues1=d.getAttributeValues();
					ArrayList<Integer> clusterIndices = d.getClusterIndex();
					
					if(i==clusterIndices.get(0))
					{
						if (d.getClassLabel() == 2){
							benignCount[i]++;
						} else{
							malignCount[i]++;
						}
						
						if(benignCount[i] > malignCount[i]){
							error[i] = (double) malignCount[i]/(benignCount[i] + malignCount[i]);
						} else{
							error[i] = (double) benignCount[i]/(benignCount[i] + malignCount[i]);

						}
						
						System.out.println(d.getscn());
						
						
					}
					
					
				}
				
				System.out.println("---------------------");
				
			}
			System.out.println(error[0]+error[1]);
			}

	
}
