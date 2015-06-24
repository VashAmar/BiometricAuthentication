import java.io.File;

import java.math.BigInteger;
import java.util.Random;
//x is client and y is server, we are trying to compute eulidean distance between the two 
//client and server values .  (x-y)^2 via computing x^2-2xy+y^2
public class vashFinal {

	static final int size = 5; //change to compute eucliean distance between of size # of elements
 	static int[] clientMsg = new int[size]; //int array of client data
	static int[] serverMsg  = new int[size]; //int array of server data
	static BigInteger[] clientCyphertextSquared = new BigInteger[size]; // biginteger array that contains x^2
	static BigInteger[] clientCyphertext = new BigInteger[size]; //biginteger array that contains x
	static BigInteger[] serverCyphertextSquared = new BigInteger[size]; //biginteger array that contains y^2
	static DamgardKey key= null;

	static final boolean debug = false;
	static final boolean getTimes = false;//setting this to true will allow you to get the times of the client and server computation but will not print euclidea distance


	public static void main(String args[])
	{
		Random r =  new Random();
		int low = 1;
		int high = 3;

		for (int i=0; i<serverMsg.length ; i++)
		{
			clientMsg[i] = r.nextInt(high-low) + low;//generate random int between high and low
			serverMsg[i] = r.nextInt(high-low) + low;

		}

		File file = new File("key.dgk");
		//key = Damgard.CreateKey(160,20,1024), message space is 2^160 key size is 1024 bits
		//Damgard.saveKeyToFile(file, key);
		key = Damgard.loadKeyFromFile(file);
		double clientTotalTime = 0;
		double serverTotalTime = 0;
		for(int i=0;i<100;i++){
		if(getTimes){
			long clientStart = System.nanoTime();
		client();
		long clientEnd = System.nanoTime(); 
		double clientTime = (clientEnd-clientStart)/1000000000.0;
		clientTotalTime += clientTime;
		long serverStart = System.nanoTime();
		server();
		long serverEnd = System.nanoTime();
		double serverTime = (serverEnd - serverStart)/1000000000.0;
		serverTotalTime += serverTime;

		}
	}
		if(!getTimes)
		{
		client();
		clientDecrypt(server());
		}

		

		s
		if(getTimes)
		{	
		System.out.println("Client Compute time  for 100 runs on " + size + " elements is: " + clientTotalTime);
		System.out.println("Server Compute time  for 100 runs on " + size + " elements is: "+ serverTotalTime);
		}
	}
	//encrypts x and x^2 and puts them in their appropriate array
	public static void client()  
	{
		int cMessage;
		int cMessageSquared;
		for(int i=0; i<serverMsg.length;i++)
		{
			cMessage = clientMsg[i]; 
			cMessageSquared = cMessage*cMessage;
			clientCyphertext[i] = Damgard.Encrypt(BigInteger.valueOf(cMessage),key);
			clientCyphertextSquared[i] = Damgard.Encrypt(BigInteger.valueOf(cMessageSquared),key);
			if(debug) System.out.println("x^2 " + i + ": " + Damgard.Decrypt(clientCyphertextSquared[i], key));
			if(debug) System.out.println("x" + i + ": " + Damgard.Decrypt(clientCyphertext[i], key));
		}
	}
	//encrypts y^2 and then calculates euclidean distance (x-y)^2
	//returns euclidean distance
	public static BigInteger server()
	{
		int sMessage;
		BigInteger tmp1 = null;
		BigInteger tmp2 = null;
		BigInteger ed = BigInteger.ONE;
		BigInteger one = BigInteger.ONE;
		BigInteger two = one.add(one);

		for(int i = 0; i<serverMsg.length; i++)
		{
			sMessage = serverMsg[i];
			sMessage = sMessage*sMessage;
			serverCyphertextSquared[i] =Damgard.Encrypt(BigInteger.valueOf(sMessage),key);
			if(debug) System.out.println("y^2 " + i + ": " + Damgard.Decrypt(serverCyphertextSquared[i], key));
		}

		for(int i = 0; i<serverMsg.length; i++)
		{
			tmp1 = clientCyphertext[i].modInverse(key.n); //tmp1  = -x
			tmp1 = tmp1.modPow(BigInteger.valueOf(serverMsg[i]),key.n); //tmp1 = -xy
			tmp1 = tmp1.modPow(two, key.n); //tmp1 = -2xy
			tmp2 = clientCyphertextSquared[i].multiply(tmp1).multiply(serverCyphertextSquared[i]); //tmp2= x^2-2xy+y^2
			tmp2 = tmp2.mod(key.n);
			ed=ed.multiply(tmp2);
			ed= ed.mod(key.n);
		}
		return ed;
		





	}
	//decrypts euclidean distance
	public static void clientDecrypt(BigInteger encrypted){
		BigInteger decrypt = encrypted;
		BigInteger dec1 = Damgard.Decrypt(decrypt, key);
		System.out.println("Euclidean Distance is: " + dec1);
		}
}





	
	


