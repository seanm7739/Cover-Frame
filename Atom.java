public class Atom implements Comparable{
	int serial;
	int neborNum;
	double x;
	double y;
	double z;
	Atom(){
		serial=0;
		neborNum=0;
		x=0;
		y=0;
		z=0;
	}
	Atom(Atom atom){
		this.serial = atom.serial;
		this.neborNum = atom.neborNum;
		this.x = atom.x;
		this.y = atom.y;
		this.z = atom.z;
	}
	public int compareTo(Object o1){
		if(this.neborNum>=((Atom)o1).neborNum){
			return 1;
	    }
		else{
			return -1;
		}
	}
}

