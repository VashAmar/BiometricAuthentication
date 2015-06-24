
import java.math.BigInteger;

public class DamgardKey implements java.io.Serializable {
	
	private static final long serialVersionUID = 1L;
	
	public final int Tbit;
	public final int Lbit;
	public final int Kbit;
	
	public final BigInteger n;
	public final BigInteger g;
	public final BigInteger h;
	public final BigInteger u;
	
	public final boolean withPrivateKey;
	public final BigInteger p;
	public final BigInteger q;
	public final BigInteger vp;
	public final BigInteger vq;
	
	public DamgardKey(int Tbit, int Lbit, int Kbit, BigInteger n, BigInteger g, BigInteger h, BigInteger u) {
		this.Tbit = Tbit;
		this.Lbit = Lbit;
		this.Kbit = Kbit;
		this.n = n;
		this.g = g;
		this.h = h;
		this.u = u;
		
		withPrivateKey = false;
		this.p = null;
		this.q = null;
		this.vp = null;
		this.vq = null;
	}
	
	public DamgardKey(int Tbit, int Lbit, int Kbit, BigInteger n, BigInteger g, BigInteger h, BigInteger u, BigInteger p, BigInteger q, BigInteger vp, BigInteger vq) {
		this.Tbit = Tbit;
		this.Lbit = Lbit;
		this.Kbit = Kbit;
		this.n = n;
		this.g = g;
		this.h = h;
		this.u = u;
		
		withPrivateKey = true;
		this.p = p;
		this.q = q;
		this.vp = vp;
		this.vq = vq;
	}
	
	@Override public String toString() {
		StringBuilder output = new StringBuilder("Damgard Key: \n");
		output.append("Tbit = "+Tbit+"\n");
		output.append("LBit = "+Lbit+"\n");
		output.append("KBit = "+Kbit+"\n");
		output.append("n = "+n+"\n");
		output.append("g = "+g+"\n");
		output.append("h = "+h+"\n");
		output.append("u = "+u+"\n");
		
		if (withPrivateKey) {
			output.append("p = "+p+"\n");
			output.append("q = "+q+"\n");
			output.append("vp = "+vp+"\n");
			output.append("vq = "+vq+"\n");
		}
		else {
			output.append("without private key");
		}

		return output.toString();
	}
}
