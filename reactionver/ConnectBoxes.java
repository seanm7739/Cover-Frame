import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ConnectBoxes {
	
	Atom atomArray[] = null;
	double atomdistanArray[][];
	double maxdistan = Double.MIN_VALUE;
	double mindiatan = Double.MAX_VALUE;
	double lr;
	//cave radius
	double cr;
	double alpha;
	double gamma;
	double rmod;
	int gbmod;
	int pushmod;
	private Scanner sca;
	
	ConnectBoxes (String pdbpath, 
				  double prl, double prh, 
				  double gx, double gy, double gz, 
				  int boxnum,
				  double lr, 
				  double beta, double alpha, int rmodule, int gbmod, double mr, int pushmod) throws IOException, InterruptedException{
		//this is multiplier of lr to nebor distances
		cr = lr * beta;
		//this is multiplier of lr to push
		this.lr = lr;
		this.alpha = alpha;
		this.gamma = 0.1;
		this.rmod = rmodule;
		this.gbmod = gbmod;
		this.pushmod = pushmod;
		//getCA from pdbqt
		getCA(pdbpath);
		getDistanceArray();
		if(rmodule == 0){
			getlrNeborNum(cr);
		}
		else if(rmodule == 1){
			getBestSDneborNum();
		}
		else if(rmodule == 2){
			getlrNeborNum(mr);
		}
		else if(rmodule == 3){
			getBestSDupbound(mr);
		}
		Arrays.sort(atomArray);
		//search good sentinel boxnum times
			System.out.println("Blue Atom:");
		for(int i = (int)(prh*(atomArray.length)) ; i < atomArray.length ; i++)
		{
			System.out.print(atomArray[i].x);
			System.out.print(" " + atomArray[i].y);
			System.out.println(" " + atomArray[i].z);
		}
			System.out.println("Surface Atom:");
		for(int i = (int)(prl*(atomArray.length)) ; i < (prh*(atomArray.length)) ; i++)
		{
			System.out.print(atomArray[i].x);
			System.out.print(" " + atomArray[i].y);
			System.out.println(" " + atomArray[i].z);
		}
			System.out.println("First Box:");
		getSentinel(1, prl, prh, gx, gy, gz);
		//System.out.println(maxdistan);
		
	}
	
	public void getCA(String path) throws IOException, InterruptedException{
		Pattern pat;
        Matcher mat;
        Boolean found;
		File f = new File(path+".pdbqt");
		sca = new Scanner(f);
		Queue<String> q = new LinkedList<String>(); 
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
	
	public void getSentinel(int boxnum, double prl, double prh , double gx, double gy, double gz){
		for(int i=0 ; i<boxnum ; i++){
			Atom ansAtom = getGoodSentinel((int)(prl*(atomArray.length)), getBlueAtom(prh, this.gbmod), gx, gy, gz);
			DecimalFormat df=new DecimalFormat("#.###");
			if(ansAtom==null){
				System.out.println("null");
				System.out.println("null");
				System.out.println("null");
				break;
			}
			else{
				System.out.print(df.format(ansAtom.x));
				System.out.print(" " + df.format(ansAtom.y));
				System.out.println(" " + df.format(ansAtom.z));
			}
		}
	}
	
	public Atom getGoodSentinel(int lindex, int hindex, double gx, double gy, double gz){
		int ans = lindex;
		//go to correct start atom
		while(ans<hindex){
			if(atomArray[ans].neborNum>0)
				break;
			ans++;
		}
		if(ans == hindex){
			return null;
		}
		//if i and j is nebor then compare i's neborNum with j
		for(int j=lindex ; j<hindex ; j++){
			double x = Math.abs(atomArray[ans].x-atomArray[j].x);
			double y = Math.abs(atomArray[ans].y-atomArray[j].y);
			double z = Math.abs(atomArray[ans].z-atomArray[j].z);
			if( x<gx && y<gy && z<gz && (atomArray[j].neborNum > atomArray[ans].neborNum) )
				ans = j;
		}
		deleteNeborNum(ans, lindex, hindex, gx, gy ,gz);
		return pushOutAns(ans, hindex);
	}
	
	public int getBlueAtom(double prh, int mod){
		if(mod == 0){
			//將範圍中，鄰居數最多者的鄰居數prh%的數量為基準，數量再其以上為藍點
			int fullnb = (int) (atomArray[atomArray.length-1].neborNum*prh);
			//System.out.println("full: "+fullnb);
			Atom t = new Atom();
			t.neborNum = fullnb;
			//System.out.println("Max: "+atomArray[atomArray.length-1].neborNum);
			int i = Arrays.binarySearch(atomArray, t);
			if(i < 0)
				i = Math.abs(i+1);
			return i;
		}
		else {
			return (int) (prh*(atomArray.length));
		}
	}
		
	public Atom pushOutAns (int ans, int hindex){
		Atom tempAtom = new Atom(atomArray[ans]);
		double[] uv= {0, 0, 0};
		//宣告blueArray，並初始化
		double[][] blueArray;
		//如果pushmod為0，代表不指定靠近的名次，則使用近gamma%的blueAtom
		if(pushmod==0){
			blueArray = new double[(int) Math.ceil(atomArray.length-hindex*gamma)][2];
		}
		//如果pushmod不為0，則使用近pushmod名的buleAtom
		else{
			blueArray = new double[pushmod][2];
		}
		for(int i=0 ; i<blueArray.length ; i++)
			blueArray[i][1]=Double.MAX_VALUE;
		//取得前gamma%或pushmod的blueAtom
		for(int i=hindex ; i<atomArray.length ; i++){
			double td = getDistance(atomArray[ans].x, atomArray[ans].y, atomArray[ans].z, atomArray[i].x, atomArray[i].y, atomArray[i].z);
			//將已有的blueAtom距離跟newd比對，較小則取代
			for(int j=0 ; j<blueArray.length ; j++){
				if(td < blueArray[j][1]){
					//紀錄atom索引，跟距離
					blueArray[j][0] = i;
					blueArray[j][1] = td;
					break;
				}
			}
		}
		//get blueAtom to ans Unit vector
		for(int i=0 ; i<blueArray.length ; i++){
			double x = atomArray[ans].x-atomArray[(int)blueArray[i][0]].x;
			double y = atomArray[ans].y-atomArray[(int)blueArray[i][0]].y;
			double z = atomArray[ans].z-atomArray[(int)blueArray[i][0]].z;
			//單位向量
			x = x/blueArray[i][1];
			y = y/blueArray[i][1];
			z = z/blueArray[i][1];
			//若pushmod為0，則加上權重
			if(pushmod==0){
				x = x/blueArray[i][1];
				y = y/blueArray[i][1];
				z = z/blueArray[i][1];
			}
			uv[0] = uv[0]+x;
			uv[1] = uv[1]+y;
			uv[2] = uv[2]+z;
		}
		
		double ud = Math.sqrt(Math.pow(uv[0], 2)+Math.pow(uv[1], 2)+Math.pow(uv[2], 2));
		//System.out.println(uv[0]+"-"+uv[1]+"-"+uv[2]);
		uv[0] = uv[0]/ud;
		uv[1] = uv[1]/ud;
		uv[2] = uv[2]/ud;
		//System.out.println(uv[0]+"-"+uv[1]+"-"+uv[2]);
		double a = lr*this.alpha;
		tempAtom.x += uv[0]*a;
		tempAtom.y += uv[1]*a;
		tempAtom.z += uv[2]*a;
		return tempAtom;
	}
	
	public void deleteNeborNum(int ans, int lindex, int hindex, double gx, double gy, double gz){
		for(int j=lindex ; j<hindex ; j++){
			double x = Math.abs(atomArray[ans].x-atomArray[j].x);
			double y = Math.abs(atomArray[ans].y-atomArray[j].y);
			double z = Math.abs(atomArray[ans].z-atomArray[j].z);
			//if i and j is nebor, delete atome
			if( x<gx && y<gy && z<gz /*&& i!=j*/){
				atomArray[j].neborNum = -1;
			}
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
	
	public void getBestSDneborNum(){
		//System.out.println(mindiatan+"-"+maxdistan);
		double SD=0;
		int continu=5; 
		double avgrArray[] = getAvgrArray(2, 1);
		for(int r=0 ; r<avgrArray.length && continu>0 ; r++){
			double avg=0;
			double sum1=0;
			double sum2=0;
			double newSD = 0;
			int tempNeborNum[] = new int[atomArray.length];
			getNeborNum(avgrArray[r], tempNeborNum);
			//compute SD
			for(int i=0 ; i<tempNeborNum.length ; i++){
				sum1 += tempNeborNum[i];
				sum2 += Math.pow(tempNeborNum[i], 2);
			}
			avg = sum1 / tempNeborNum.length;
			newSD = Math.sqrt( (sum2/tempNeborNum.length)-Math.pow(avg, 2) );
			//get better neborNum if SD is better
			if (newSD > SD){
				SD = newSD;
				continu = 5;
				for(int i=0 ; i<atomArray.length ; i++)
					atomArray[i].neborNum = tempNeborNum[i];
				this.cr = avgrArray[r];
			}
			else{
				continu--;
			}	
		}
	}

        public void getBestSDupbound(double mr){
                //System.out.println(mindiatan+"-"+maxdistan);
                double SD=0;
                int continu=5;
                double avgrArray[] = getrArray(mr);
                for(int r=0 ; r<avgrArray.length && continu>0 ; r++){
                        double avg=0;
                        double sum1=0;
                        double sum2=0;
                        double newSD = 0;
                        int tempNeborNum[] = new int[atomArray.length];
                        getNeborNum(avgrArray[r], tempNeborNum);
                        //compute SD
                        for(int i=0 ; i<tempNeborNum.length ; i++){
                                sum1 += tempNeborNum[i];
                                sum2 += Math.pow(tempNeborNum[i], 2);
                        }
                        avg = sum1 / tempNeborNum.length;
                        newSD = Math.sqrt( (sum2/tempNeborNum.length)-Math.pow(avg, 2) );
			//System.out.println(" " + newSD + " " + avgrArray[r] + " " + mr + " " + continu);
                        //get better neborNum if SD is better
                        if (newSD > SD){
                                SD = newSD;
				continu = 5;
                                for(int i=0 ; i<atomArray.length ; i++)
                                        atomArray[i].neborNum = tempNeborNum[i];
                                this.cr = avgrArray[r];
				if(this.cr > mr){
					//this.cr = mr;
					getlrNeborNum(mr);
					continu=0;
				}
                        }
                        else{
                                continu--;
                        }
                }
        }

        private double[] getrArray(double kh){
                int h = (int) (kh+6);
                double avgrArray[] = new double[h*2];
                for(int t=0; t<h*2 ; t++){
                	avgrArray[t] = t*0.5 + 3;
			//System.out.println(t);
                }
                return avgrArray;
        }
	
	private double[] getAvgrArray(double kl, double kh){
		int l = (int) (atomdistanArray[0].length/kl);
		int h = (int) (atomdistanArray[0].length/kh);
		double avgrArray[] = new double[h-l];
		for(int j=l, t=0; j<h ; j++, t++){
			for(int i=0 ; i<atomdistanArray.length ; i++){
				avgrArray[t]+=atomdistanArray[i][j];
			}
			avgrArray[t] = avgrArray[t] / (atomdistanArray.length-1);
		}
		return avgrArray;
	}
	
	public void getlrNeborNum(double cr){
		int tempNeborNum[] = new int[atomArray.length];
		//getNeborNum
		getNeborNum(cr, tempNeborNum);
		//output NeborNum in atomArray.neborNum;
		for(int i=0 ; i<atomArray.length ; i++)
			atomArray[i].neborNum = tempNeborNum[i];
	}
	
	private void getNeborNum(double cr, int tempArray[]){
		for(int i=0 ; i<tempArray.length ; i++){
			int j = Arrays.binarySearch(atomdistanArray[i], cr);
			j = Math.abs(j+1)-1;
			tempArray[i] = j;	
		}
	}
	
	private double getDistance (double x1, double y1, double z1, double x2, double y2, double z2){
		return Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2)+Math.pow(z2-z1, 2));
	}

	//[0]=pdb, [1]=prl, [2]=prh, [3,4,5]=gx,y,z, [6]=box number, [7]=ligen radius, [8]=beta, [9]=alpha
	//[10]=rmodule, [11]=gbmod, [12]=mr, [13]=pushmode
	public static void main(String[] args) throws IOException, NumberFormatException, InterruptedException{
		new ConnectBoxes(args[0], 
						 Float.parseFloat(args[1]), 
						 Float.parseFloat(args[2]),
						 (Float.parseFloat(args[3])/2), 
						 (Float.parseFloat(args[4])/2), 
						 (Float.parseFloat(args[5])/2), 
						 Integer.parseInt(args[6]),
						 Float.parseFloat(args[7]),
						 Float.parseFloat(args[8]),
						 Float.parseFloat(args[9]), 
						 Integer.parseInt(args[10]),	//rmodule
						 Integer.parseInt(args[11]),	//get blue atom mod
						 Float.parseFloat(args[12]),	//manual radius
						 Integer.parseInt(args[13])		//push mod
						);
	}
}


