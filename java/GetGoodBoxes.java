import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.text.html.HTMLDocument.HTMLReader.HiddenAction;

public class GetGoodBoxes {
	
	Atom atomArray[] = null;
	double atomdistanArray[][];
	double maxdistan = Double.MIN_VALUE;
	double mindiatan = Double.MAX_VALUE;
	
	public void showArray(double[] t){
		for(double i:t){
			System.out.println(i+"\t");
		}
	}
	public void showArray(double[][] t){
		for(double[] i:t){
			for(double j:i){
				System.out.print(j+"\t");
			}
			System.out.println();
		}
	}

	
	GetGoodBoxes (String pdbpath, double lr, double percentl, double percenth, double gx, double gy, double gz, int boxnum) throws IOException, InterruptedException{
		//getCA form pdbqt by perl
		getCA(pdbpath);

		getDistanceArray();
		//showArray(atomdistanArray);
		getNeborNumResult(lr);
	
		Arrays.sort(atomArray);
		
		//search good sentinel boxnum times
		getSentinel(boxnum, percentl, percenth, gx, gy, gz, lr);
		/*
		for(int i=0 ; i<atomArray.length ; i++){
			System.out.println(atomArray[i].serial+"\t"+atomArray[i].neborNum+"\t"+atomArray[i].x+"\t"+atomArray[i].y+"\t"+atomArray[i].z);
		}*/
		
	}
	
	public void getSentinel(int boxnum, double percentl, double percenth , double gx, double gy, double gz, double lr){
		for(int i=0 ; i<boxnum ; i++){
			int ans = getGoodSentinel((int)(percentl*(atomArray.length)), (int)(percenth*(atomArray.length)), gx, gy, gz ,lr);
			//String a = atomArray[ans].x+"\t"+atomArray[ans].y+"\t"+atomArray[ans].z;
			System.out.println(atomArray[ans].x);
			System.out.println(atomArray[ans].y);
			System.out.println(atomArray[ans].z);
		}
	}
	
	public int getGoodSentinel(int lindex, int hindex, double gx, double gy, double gz, double lr){
		//if i and j is nebor then compare i's neborNum with j
		/*for(int t=i ; t<hindex ; t++){
			System.out.println(atomArray[t].serial+"\t"+atomArray[t].neborNum+"\t"+atomArray[t].x+"\t"+atomArray[t].y+"\t"+atomArray[t].z);
		}*/
		int ans = lindex;
		//go to correct start atom
		while(ans<hindex){
			if(atomArray[ans].neborNum>0 || ans==(atomArray.length-1)){
				break;
			}
			ans++;
		}
		
		for(int j=lindex ; j<hindex ; j++){
			double x = Math.abs(atomArray[ans].x-atomArray[j].x);
			double y = Math.abs(atomArray[ans].y-atomArray[j].y);
			double z = Math.abs(atomArray[ans].z-atomArray[j].z);
			if( x<gx && y<gy && z<gz && (atomArray[j].neborNum > atomArray[ans].neborNum) )
				ans = j;
		}
		pushOutAns(ans, hindex, lr);
		deleteNeborNum(ans, lindex, hindex, gx, gy ,gz);
		return ans;
	}
	
	public void pushOutAns (int ans, int hindex, double lr){
		double[] uv= {0, 0, 0};
		double x1 = atomArray[ans].x;
		double y1 = atomArray[ans].y;
		double z1 = atomArray[ans].z;
		//get BlueAtom nebor by lr
		for(int i=hindex ; i<atomArray.length ; i++){
			double td = getDistance(x1, y1, z1, atomArray[i].x, atomArray[i].y, atomArray[i].z);
			if(td < lr){
				//get all blueNebor to ans Unit vector
				double x = x1-atomArray[i].x;
				double y = y1-atomArray[i].y;
				double z = z1-atomArray[i].z;
				x = x/td;
				y = y/td;
				z = z/td;
				uv[0] = uv[0]+x;
				uv[1] = uv[1]+y;
				uv[2] = uv[2]+z;
			}
		}
		double td = Math.sqrt(Math.pow(uv[0], 2)+Math.pow(uv[1], 2)+Math.pow(uv[2], 2));
		uv[0] = uv[0]/td;
		uv[1] = uv[1]/td;
		uv[2] = uv[2]/td;
		//get all Unit vector sum and push ans by n
		atomArray[ans].x += uv[0]*2;
		atomArray[ans].y += uv[1]*2;
		atomArray[ans].z += uv[2]*2;
	}
	
	public void deleteNeborNum(int i, int lindex, int hindex, double gx, double gy, double gz){
		for(int j=lindex ; j<hindex ; j++){
			double x = Math.abs(atomArray[i].x-atomArray[j].x);
			double y = Math.abs(atomArray[i].y-atomArray[j].y);
			double z = Math.abs(atomArray[i].z-atomArray[j].z);
			//if i and j is nebor, delete atome
			if( x<gx && y<gy && z<gz /*&& i!=j*/){
				atomArray[j].neborNum = -1;
			}
		}
	}
	
	
	public void getCA(String path) throws IOException, InterruptedException{
		Pattern pat;
        Matcher mat;
        Boolean found;
		File f = new File(path+".pdbqt");
		Scanner sca = new Scanner(f);
		Queue q = new LinkedList(); 
	    pat = Pattern.compile("ATOM +[0-9]+ +CA");
		while(sca.hasNext()){
			String t = sca.nextLine();
			mat = pat.matcher(t);
			found = mat.find();
			if(found){
				q.offer(t);
			}
		}
	    atomArray = new Atom[q.size()];
	    for(int i=0 ; i<atomArray.length ; i++){
        	atomArray[i] = new Atom();
        	String t[] = ((String)q.poll()).split(" +");
        	atomArray[i].serial = Integer.parseInt(t[1]);
        	atomArray[i].x = Double.parseDouble(t[6]);
        	atomArray[i].y = Double.parseDouble(t[7]);
        	atomArray[i].z = Double.parseDouble(t[8]);
        }
	}
	
	public void getDistanceArray(){
		atomdistanArray = new double[atomArray.length][atomArray.length];
		for(int i=0 ; i<atomArray.length ; i++){
			for (int j=i ; j<atomArray.length ; j++){
				atomdistanArray[i][j] = getDistance(atomArray[i].x, atomArray[i].y, atomArray[i].z, 
			               atomArray[j].x, atomArray[j].y, atomArray[j].z);
				atomdistanArray[j][i] = atomdistanArray[i][j];
				if(atomdistanArray[i][j] > maxdistan)
					maxdistan = atomdistanArray[i][j];
				if(atomdistanArray[i][j] < mindiatan && atomdistanArray[i][j] != 0)
					mindiatan = atomdistanArray[i][j];
			}
			Arrays.sort(atomdistanArray[i]);
		}
	}
	
	public void getNeborNumResult(double lr){
		//System.out.println(mindiatan+"-"+maxdistan);
		/*
		double SD=0;
		double avg=0;
		double sum1=0;
		double sum2=0;
		*/
		int tempNeborNum[] = new int[atomArray.length];
		//getNeborNum
		getNeborNum(lr, tempNeborNum);
		/*
		//getSD
		for(int i=0 ; i<tempNeborNum.length ; i++){
			sum1 += tempNeborNum[i];
			sum2 += Math.pow(tempNeborNum[i], 2);
		}
		avg = sum1 / tempNeborNum.length;
		SD = Math.sqrt( (sum2/tempNeborNum.length)-Math.pow(avg, 2) );
		//System.out.println(lr+", "+SD);
		*/
		//output NeborNum in atomArray.neborNum;
		for(int i=0 ; i<atomArray.length ; i++)
			atomArray[i].neborNum = tempNeborNum[i];
	}
	
	private void getNeborNum(double lr, int tempArray[]){
		for(int i=0 ; i<tempArray.length ; i++){
			int j = Arrays.binarySearch(atomdistanArray[i], lr);
			j = Math.abs(j+1)-1;
			tempArray[i] = j;	
		}
	}
	
	private double getDistance (double x1, double y1, double z1, double x2, double y2, double z2){
		return Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2)+Math.pow(z2-z1, 2));
	}
	

	
	public static void main(String[] args) throws IOException, NumberFormatException, InterruptedException{
		                //	
		new GetGoodBoxes(args[0], Float.parseFloat(args[1]), Float.parseFloat(args[2]), Float.parseFloat(args[3]),
				(Float.parseFloat(args[4])/2), (Float.parseFloat(args[5])/2), (Float.parseFloat(args[6])/2), 
				Integer.parseInt(args[7]));
	}
}


