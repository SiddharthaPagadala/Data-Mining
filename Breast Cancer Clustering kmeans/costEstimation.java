package BreastCancerKMeans;

import java.util.*;

public class costEstimation {

	public static void main(String[] args) {
		int biopsysCost = 0, masectomyCost=0;
		Random r =new Random();
		
		for(int i=1; i<700;i++){
			biopsysCost = biopsysCost + r.nextInt((5000-1000)+1) + 1000;
		}
		
		for(int i=1;i<242;i++){
			masectomyCost = masectomyCost+r.nextInt((55000-15000)+1) + 15000;
		}
		System.out.println("Total biopsy cost:"+biopsysCost);
		System.out.print("Total masectomy cost:"+masectomyCost);
	}

}
