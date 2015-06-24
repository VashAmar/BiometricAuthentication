
import java.math.BigInteger;
import java.util.Arrays;

public class DForComparisonResult
{
    public BigInteger d;
    public BigInteger[] d_i;
    
    @Override public String toString() {
    	StringBuilder builder = new StringBuilder();
    	builder.append("d: "+d+"\n");
    	builder.append("d_i: "+Arrays.toString(d_i));
    	
    	return builder.toString();
    }
}