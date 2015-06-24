
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OptionalDataException;
import java.io.OutputStream;
import java.io.StreamCorruptedException;
import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

//import android.util.Log;

public class Damgard {
	/*
	 * Lbit is the plaintext length in bits.
	 * 
	 * Tbit, Kbit are security parameters. Length of RSA modulus is 2*Kbit.
	 * http://www.keylength.com/en/4/
	 * security parameter combinations:
	 * T: 160	224
	 * K: 512	1024
	 */
	
	public static final BigInteger minusOne = BigInteger.ONE.negate();
	public static final BigInteger zero = BigInteger.ZERO;
	public static final BigInteger one = BigInteger.ONE;
	public static final BigInteger two = one.add(one);
	public static final BigInteger three = two.add(one);
	
	public static Map<Integer, BigInteger[]> powersOf2 = new HashMap<Integer, BigInteger[]>();
	public static Map<DamgardKey, BigInteger[]> gvpps = new HashMap<DamgardKey, BigInteger[]>();
	
	public static long timeForAlice = 0;
	public static long timeForBob = 0;
	
	public static DamgardKey CreateKey(int Tbit, int Lbit, int Kbit) {
		Random random = new SecureRandom();
		
		//choose vp, vq as primes of length Tbit
		BigInteger vp = BigInteger.probablePrime(Tbit, random);
		BigInteger vq = BigInteger.probablePrime(Tbit, random);
		
		//set u to 2 ^ (Lbit)
		BigInteger u = BigInteger.valueOf(2).pow(Lbit);
		
		// create p as a prime equal to (RandomPrime * u * vp) + 1
		BigInteger p;
		BigInteger f1;
		int count=0;
	   	while (true) {
			count++; 
			f1 = BigInteger.probablePrime(Kbit-Tbit-Lbit+1, random);	 
			p = vp.multiply(u).multiply(f1);
			p = p.add(one); 
			if (p.isProbablePrime(10)) {
				break;
			}
	   	}
	   	System.out.println("number of trials for p: "+count);
	   	
	   	// create q as a prime equal to RandomPrime * u * vq
		BigInteger q;
		BigInteger f2;
		count=0;
	   	while (true) {
			count++; 
			f2 = BigInteger.probablePrime(Kbit-Tbit-Lbit+1, random);	 
			q = vq.multiply(u).multiply(f2);
			q = q.add(one); 
			if (q.isProbablePrime(10)) {
				break;
			}
	   	}
	   	System.out.println("number of trials for q: "+ count);	   	
	   	
	    //compute n as the product of p * q
		BigInteger n = p.multiply(q);
		
		// choose a number from Z^*_p and check its order
		//exp1=(2^(Lbit - 1) * vp * f1)
		//exp2=u * vp
		//exp3=u *f1
		BigInteger xp;
		count=0;
		int num_bits;
		BigInteger exp1 = BigInteger.valueOf(2).pow(Lbit-1).multiply(vp).multiply(f1); 
		BigInteger exp2 = u.multiply(vp);
		BigInteger exp3 = u.multiply(f1);
		while(true) {
			count++;
			num_bits = p.bitLength();
			xp = new BigInteger(num_bits, random); //TODO: CHECK xp BITSIZE!
			if (xp.compareTo(p) < 0) {
				//tmp=xp ^exp1 mod p
				BigInteger tmp = xp.modPow(exp1, p);
				if (tmp.compareTo(one)==0) {
					continue;
				}

				//tmp=xp^exp2 mod p
				tmp = xp.modPow(exp2, p);
				if (tmp.compareTo(one)==0) {
					continue;
				}

				//tmp=xp ^ exp3 mod p
				tmp = xp.modPow(exp3, p);
				if (tmp.compareTo(one)==0) {
					continue;
				}

				break;
			}
		}
		System.out.println("number of trials for xp: "+ count);
		
		// choose a number from Z^*_q and check its order
		//exp1=(2^(Lbit - 1) * vq * f2)
		//exp2=u * vq
		//exp3=u *f2
		BigInteger xq;
		count = 0;
	   	exp1 = BigInteger.valueOf(2).pow(Lbit - 1).multiply(vq).multiply(f2);
		exp2 = u.multiply(vq); 
		exp3 = u.multiply(f2);    	
		while(true){
			count++;
		    num_bits=q.bitLength();
		    xq = new BigInteger(num_bits, random);	
			//random_mpz(xq, Kbit/8); 
			if (xq.compareTo(q) < 0){
				//tmp=xp ^exp1 mod q
				BigInteger tmp = xq.modPow(exp1, q);
				if (tmp.compareTo(one) == 0){
					continue; 
				}	
				//tmp=xp^exp2 mod q
				tmp = xq.modPow(exp2, q);
				if (tmp.compareTo(one) == 0){
					continue; 
				}
				//tmp=xp ^ exp3 mod q
				tmp = xq.modPow(exp3, q); 
				if (tmp.compareTo(one) == 0){
					continue; 
				}

				break; 
			}
		}
		System.out.println("number of trials for xq: "+ count);

		//compute g
		//first, compute CRT: g = xp*q*(q^{-1} mod p) + xq*p*(p^{-1} mod q) mod n
		BigInteger firstSummand = q.modInverse(p).multiply(q).mod(n).multiply(xp).mod(n);
		BigInteger secondSummand = p.modInverse(q).multiply(p).mod(n).multiply(xq).mod(n);
		BigInteger g = firstSummand.add(secondSummand).mod(n);
		
		//second, compute g = g ^(f1f2) mod n...so it will be of uvpvq order...
		BigInteger exponent = f1.multiply(f2); 
		g = g.modPow(exponent, n); 
		
		//compute h 
		BigInteger h;
		while(true){
			num_bits=n.bitLength();
			h = new BigInteger(num_bits, random);
			if (h.compareTo(n) < 0){
				break; 
			}
		}
		exponent = f1.multiply(f2).multiply(u);	
		h = h.modPow(exponent, n);
		
		DamgardKey key = new DamgardKey(Tbit, Lbit, Kbit, n, g, h, u, p, q, vp, vq);
		
		// prepare pow2s and gvpps for decryption
		getPowersOf2(Lbit);
		getGvpp(key);
		
		return key;
	}
	
	public static BigInteger Encrypt(BigInteger m, DamgardKey key){	
		Random random = new SecureRandom();
		BigInteger r = new BigInteger((int) (2.5*key.Tbit), random);
		
		// EncryptedMessage = h^r * g^m (mod n)
		BigInteger first = key.h.modPow(r, key.n);
		BigInteger second = key.g.modPow(m, key.n);
		
		return first.multiply(second).mod(key.n); 
	}	
	
	public static BigInteger EncryptDouble(double m, DamgardKey key, int precision) {
		double shift = Math.pow(10, precision);
		long message = Math.round(m*shift);
		return Encrypt(BigInteger.valueOf(message), key);
	}
	
	public static double DecryptDouble(BigInteger em, DamgardKey key, int precision, boolean isProduct) {
		BigInteger message = Damgard.Decrypt(em, key);
		//System.out.println(message.toString());
		
		if (isProduct) {
			BigInteger shift = BigInteger.TEN.pow(precision);
			message = message.divide(shift);
		}
		
		return message.longValue()/Math.pow(10, precision);		
	}
	
	/**
	 * homomorphic multiplication on encrypted double value
	 * @return
	 */
	public static BigInteger hmul(BigInteger em, double numberB, DamgardKey key, int precision) {
		BigInteger exponent = BigInteger.valueOf(Math.round(numberB*Math.pow(10, precision)));
		return em.modPow(exponent, key.n);
	}
	
	public static BigInteger Decrypt(BigInteger em, DamgardKey key) {
		if (! key.withPrivateKey) {
			throw new IllegalArgumentException("private key missing");
		}
		
		int i; 
		int xi[] = new int[key.Lbit]; 
		BigInteger pow2[] = getPowersOf2(key.Lbit);
		BigInteger gvpp[] = getGvpp(key);
	
		BigInteger y;
		BigInteger yi;
		
		y = em.modPow(key.vp, key.p);

		for(i=0; i < key.Lbit; i ++) {
			yi = y.modPow(pow2[key.Lbit-1-i], key.p);
			if (yi.compareTo(one) == 0) {
				xi[i]=0;
			}
			else{
				xi[i]=1;
				y = y.multiply(gvpp[i]).mod(key.p);
			}		
		}
	
		BigInteger m = BigInteger.valueOf(xi[0]);
		for(i=1; i < key.Lbit; i++) {
			if (xi[i] == 1) {
				m = m.add(pow2[i]);
			}
		}	
		
		return m;
	}
	
	/**
	 * Returns true if em is an encryption of zero under the key. 
	 * This method is much faster than Decrypt.
	 *  
	 * @param e_value
	 * @param key
	 * @return
	 */
	public static boolean IsEncryptedZero(BigInteger em, DamgardKey key) {
		// tmp = em^(vp*vq) mod n
		 BigInteger tmp = em.modPow(key.vp, key.n).modPow(key.vq, key.n);
		
		if (tmp.equals(one)) {
			return true; 
		}
		
		return false;
	}
	
	// return encrypted*(h^r) (mod n)
	public static BigInteger ReRandomize(BigInteger encrypted, DamgardKey key) {
		Random random = new SecureRandom();
		
		BigInteger exponent = new BigInteger((int) (2.5*key.Tbit), random);
		BigInteger randomizer = key.h.modPow(exponent, key.n);
		
		return encrypted.multiply(randomizer).mod(key.n);
	}
	
	/*
	 * Based on Erkin et. al.'s Privacy-preserving face recognition, Section 5
	 *
	 * Compare the size of encrypted values a and b. 
	 * M is an upper bound for a, b sizes in bits.
	 * kappa is a security parameter.
	 * (M+kappa) must be smaller than key.Lbit! 
	 * @return -1 if a is smaller than b, +1 otherwise
	 */
	public static int CompareEncrypted(BigInteger a, BigInteger b, BigInteger M, BigInteger kappa, DamgardKey key) {
		Random random = new SecureRandom();
		
		StartCounting(); 
		long start = System.currentTimeMillis();
		
		// z = 2^M + a - b
		BigInteger twoToM = BigInteger.valueOf(2).pow(M.intValue());
		BigInteger e_twoToM = Damgard.Encrypt(twoToM, key);
		BigInteger e_z = e_twoToM.multiply(a).mod(key.n).multiply(b.modInverse(key.n)).mod(key.n);
		
		// d = z + r
		BigInteger r = new BigInteger(M.add(kappa).intValue(), random);
		BigInteger e_r = Encrypt(r, key);
		BigInteger e_d = e_z.multiply(e_r).mod(key.n);
		e_d = ReRandomize(e_d, key);

		timeForBob +=System.currentTimeMillis()-start;
		start = System.currentTimeMillis();
		
		// reduce d to (d mod 2^M)
		DForComparisonResult dForComparison = CompareEncrypted_GetDMod2M_DBits_Encrypted(e_d, M, key);
		
		timeForAlice += System.currentTimeMillis()-start;
		start = System.currentTimeMillis();

		BigInteger e_d_hat = dForComparison.d;
		BigInteger[] e_d_bits = dForComparison.d_i;
		
		// z_tilde = (d mod 2^M) - (r mod 2^M)
		BigInteger r_hat = r.mod(twoToM);
		BigInteger e_r_hat = Encrypt(r_hat, key);
		BigInteger e_z_tilde = e_d_hat.multiply(e_r_hat.modInverse(key.n)).mod(key.n);
		
		// TODO: remove debug
//		System.out.println("twoToM: "+twoToM);				
//		System.out.println("e_twoToM decrypted: "+Damgard.Decrypt(e_twoToM, key));				
//		System.out.println("e_z decrypted: "+Damgard.Decrypt(e_z, key));	
//		System.out.println("e_d decrypted: "+Damgard.Decrypt(e_d, key));
//		System.out.println("e_d_hat decrypted: "+Damgard.Decrypt(e_d_hat, key));	
//		System.out.println("e_r decrypted: "+Damgard.Decrypt(e_r, key));
//		System.out.println("e_r_hat decrypted: "+Damgard.Decrypt(e_r_hat, key));				
//		System.out.println("e_z_tilde decrypted: "+Damgard.Decrypt(e_z_tilde, key));		
		
		// if (r mod 2^M) < (d mod 2^M), z_tilde is the value of (z mod 2^M)
		// we have to find out if it is the case.
		
		// for technical reasons:  r̂ = 2*r̂, d_hat = 2*d_hat+1
		// modified d's bits are already in e_d_bits
		BigInteger r_hat_comparison = r_hat.multiply(two);
		boolean[] r_bits = new boolean[M.intValue()+1];
		for (int i = 0; i < M.intValue()+1; i++) {
			r_bits[i] = r_hat_comparison.testBit(i) ? true : false;
		}
		
		timeForBob+=System.currentTimeMillis()-start;
		BigInteger e_lambda = CompareBitArrays(e_d_bits, r_bits, key);
		start = System.currentTimeMillis();
		
		// fix possible underflow
		// (z mod 2^M) = z_tilde + lambda*2^M 
		BigInteger e_z_mod = e_z_tilde.multiply(e_lambda.modPow(twoToM, key.n)).mod(key.n);
		
		// z[l] = 2^(-M)*(z - (z mod 2^M))
		// ! it is ok to just use (z - (z mod 2^M)), we get 2^M*(z_l)
		BigInteger e_z_l = e_z.multiply(e_z_mod.modInverse(key.n));

		// TODO: remove debug
		//System.out.println("e_z_mod decrypted: "+Damgard.Decrypt(e_z_mod, key));			
		//System.out.println("e_z_l decrypted: "+Damgard.Decrypt(e_z_l, key));
		
		timeForBob += System.currentTimeMillis()-start;
		start = System.currentTimeMillis();
		
		// Alice can decrypt the result and find out, if a<b or not.
		boolean isEncryptedZero = IsEncryptedZero(e_z_l, key);
		
		timeForAlice += System.currentTimeMillis()-start;
		
		if (isEncryptedZero) {
			return -1;
		} else {
			return 1;
		}
	}
	
	private static void StartCounting() {
		timeForAlice = 0;
		timeForBob = 0;
	}
	
	/**
	 * Based on Erkin et. al.'s Privacy-preserving face recognition, Section 5
	 * decrypt d, reduce it modulo 2^M, encrypt it again
	 * 
	 * DO NOT CALL DIRECTLY! Does ALICE's job for SubSection 5.2
	 */
	public static DForComparisonResult CompareEncrypted_GetDMod2M_DBits_Encrypted(BigInteger d, BigInteger M, DamgardKey key) {
		BigInteger decrypted = Damgard.Decrypt(d, key);
		BigInteger dMod2M = decrypted.mod(two.pow(M.intValue()));
		
		BigInteger dModForComparison = dMod2M.multiply(two).add(one);
		
		BigInteger[] d_i = new BigInteger[M.intValue()+1];
		for (int i = 0; i < M.intValue()+1; i++) {
			boolean bitValue = dModForComparison.testBit(i) ? true : false;
			d_i[i] = Damgard.Encrypt(bitValue ? one : zero, key);
		}
		
		DForComparisonResult dWithBitsEncrypted = new DForComparisonResult();
		dWithBitsEncrypted.d = Damgard.Encrypt(dMod2M, key);
		dWithBitsEncrypted.d_i = d_i;
		
		return dWithBitsEncrypted;
	}
	
	public static boolean[] GetBitArray(BigInteger a, int M) {
		boolean[] result = new boolean[M];
		for (int i = 0; i < M; i++) {
			result[i] = a.testBit(i);
		}
		return result;
	}
	
	/**
	 * Based on SubSection 5.3 from Etkin et al's Privacy-Preserving face recognition
	 * 
	 * Compare two BigIntegers d, r represented as bit arrays of the same length. 
  	 * First of the two arrays is encrypted using key.
  	 * 
  	 * The result is:
  	 *   an encryption of 1 if r > d
  	 *   an encryption of 0 otherwise 
	 */
	public static BigInteger CompareBitArrays(BigInteger[] e_d_bits, boolean[] r_bits, DamgardKey key) {
		if (e_d_bits.length != r_bits.length) 
			throw new IllegalArgumentException("bit array lengths for d and r do not match!");
		
		long start = System.currentTimeMillis();
		
		// encrypt r_bits and minus r_bits
		BigInteger[] e_r_bits = new BigInteger[r_bits.length];
		BigInteger[] e_mr_bits = new BigInteger[r_bits.length];
		for (int i = 0; i < r_bits.length; i++) {
			BigInteger bitValue = r_bits[i] ? one : zero;
			BigInteger minusBitValue = r_bits[i] ? minusOne : zero;
			e_r_bits[i] = Encrypt(bitValue, key);
			e_mr_bits[i] = Encrypt(minusBitValue, key);
		}

		// w[i] = d[i] XOR r[i]
		BigInteger[] e_w = new BigInteger[r_bits.length];
		for (int i = 0; i < r_bits.length; i++) {
			if (r_bits[i]==true) {
				BigInteger e_one = Encrypt(one, key);
				e_w[i] = e_one.multiply(e_d_bits[i].modInverse(key.n)).mod(key.n);
			} else {
				e_w[i] = e_d_bits[i];
			}
		}
		
		// sumw[i] = Sum_{j=i+1}^{M}(w[i])
		BigInteger[] e_sumw = new BigInteger[e_w.length];
		e_sumw[e_w.length-1] = Encrypt(zero, key);
		for (int i = e_d_bits.length-2; i >= 0; i--) {
			e_sumw[i] = e_sumw[i+1].multiply(e_w[i+1]).mod(key.n);
		}
		
		// choose s \in {-1, +1} randomly
		Random random = new SecureRandom();
		BigInteger s = random.nextBoolean() ? minusOne : one;
		BigInteger e_s = Encrypt(s, key);
		
		//TODO: remove debug
		//System.out.println("e_s decrypted: "+Damgard.Decrypt(e_s, key));				
		
		// c[i]  = d_i - r_i + s + 3*sumw[i]
		BigInteger[] e_c = new BigInteger[r_bits.length];
		for (int i = 0; i < r_bits.length; i++) {
			e_c[i] = e_d_bits[i].multiply(e_mr_bits[i]).mod(key.n).multiply(e_s).mod(key.n).multiply(e_sumw[i].modPow(three, key.n)).mod(key.n);
		}
		
		// mask each e_c[i] with random even t[i] (that is, t[i] is coprime to 2^Lbit) 
		// and rerandomize e_c[i]
		BigInteger[] t = new BigInteger[r_bits.length];
		for (int i = 0; i < r_bits.length; i++) {
			t[i] = new BigInteger(key.Lbit-1, random).setBit(0);
			e_c[i] = e_c[i].modPow(t[i], key.n);
			e_c[i] = ReRandomize(e_c[i], key);
		}
		
		timeForBob += System.currentTimeMillis()-start;
		start = System.currentTimeMillis();
		
		// ask Alice, if e_c contains encrypted zero. answer is encrypted
		BigInteger e_lambda_hat = DoesContainEncryptedZero(e_c, key);
		BigInteger e_lambda;	
		
		timeForAlice += System.currentTimeMillis()-start;
		start = System.currentTimeMillis();
		
		// lambda = lambda_hat XOR (!s)
		// that is, lambda = 1-lambda_hat if s==-1
		if (s.equals(one)) {
			e_lambda = e_lambda_hat;
		} else {
			e_lambda = Encrypt(one, key).multiply(e_lambda_hat.modInverse(key.n));
		}

		//TODO: remove debug
		//System.out.println("e_lambda_hat decrypted: "+Damgard.Decrypt(e_lambda_hat, key));		
		//System.out.println("e_lambda decrypted: "+Damgard.Decrypt(e_lambda, key));				
		
		timeForBob += System.currentTimeMillis() - start;
		return e_lambda;
	}
	
	
	/*
	 * Given an array of encrypted values, return:
	 *    Encrypt(one) if one of those values contains zero.
	 *    Encrypt(zero) if all values are   
	 */
	public static BigInteger DoesContainEncryptedZero(BigInteger[] e_array, DamgardKey key) {
		for (BigInteger e_value : e_array) {
			if (IsEncryptedZero(e_value, key)) {
				return Encrypt(one, key);
			}
		}
		return Encrypt(zero, key);
	}


	public static BigInteger[] getPowersOf2(int L) {
		if (powersOf2.get(L)!=null) {
			return powersOf2.get(L);
		}
		
		BigInteger[] pow2 = new BigInteger[L];
		pow2[0]=one;
		pow2[1]=one.add(one);
		for (int i = 2; i<L; i++) {
			pow2[i] = pow2[i-1].multiply(pow2[1]);
		}
		powersOf2.put(L, pow2);
		return pow2;
	}
	
	public static BigInteger[] getGvpp(DamgardKey key) {
		if (gvpps.get(key)!=null) {
			return gvpps.get(key);
		}
		
		BigInteger[] gvpp = new BigInteger[key.Lbit];
		BigInteger[] pow2 = getPowersOf2(key.Lbit);
		
		BigInteger gvp = key.g.modPow(key.vp, key.p);
		BigInteger tmp = key.u.subtract(one);
		
		for (int i = 0; i<key.Lbit; i++) {
			gvpp[i] = gvp.modPow(pow2[i], key.p);
			gvpp[i] = gvpp[i].modPow(tmp, key.p);
		}
		gvpps.put(key, gvpp);
		return gvpp;
	}
	
	public static DamgardKey loadKeyFromFile(File file) {
		if (! file.exists()) {
			throw new IllegalArgumentException("file "+file+" does not exist");
		}
		
		try {
			ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream(file));
			DamgardKey key = (DamgardKey) inputStream.readObject();
			Damgard.getGvpp(key);
			Damgard.getPowersOf2(key.Lbit);
			
			return key;
		} catch (Exception e) {
			//System.out.println("LOADKEY","load key from file error"+file, e);
			throw new RuntimeException("key load error");
		}
	}
	
	public static void saveKeyToFile(File file, DamgardKey key) {
		try {
			OutputStream fos = new FileOutputStream(file);
			ObjectOutputStream outputStream = new ObjectOutputStream(fos);
			outputStream.writeObject(key);
			fos.close();
		} catch (Exception e) {
			//System.out.println("SAVEKEY","save key fo file error "+file, e);
			throw new RuntimeException("save key error", e);
		}
	}
}