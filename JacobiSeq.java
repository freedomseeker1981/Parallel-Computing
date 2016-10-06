
/*
 * Name : Reda Elbahi
 * Course : Parallel computing I
 * Prof: Alan Kaminsky
 * Project: Programming Project I sequential version
 * Description: This program implements the jacobi method to solve linear algebra equations and test the convergnces
 * against the value of 10^-8 and stops the iteration at that level and prints the results
 */
import edu.rit.util.Random;
import edu.rit.pj.Comm;

public class JacobiSeq {


	public static void main(String args[]) throws Exception{

		
		long t1 = System.currentTimeMillis();
		
		if(args.length != 2)
		{
			System.err.println("worng input :   Please check your argument ");

		}
		
		int n = 0;
		double[] y = null;
		double[] x = null;
		double[][] A ={{8.0,2.0,3.0},{-1.0 ,6.0 ,4.0},{5.0 ,1.0, 7.0}} ;
		
		double[] b = null;

		try{
			long seed = Long.parseLong(args[1]);
			n = Integer.parseInt(args[0]);
			A = new double [n] [n];
			b = new double [n];
			Random prng = Random.getInstance (seed);
			for (int i = 0; i < n; ++ i)
			{
				for (int j = 0; j < n; ++ j)
					A[i][j] = prng.nextDouble()*9.0 + 1.0;
				A[i][i] += 10.0*n;
				b[i] = prng.nextDouble()*9.0 + 1.0;
			}

			x = new double [n];//pre-condition
			y = new double [n];// post-condition
			// set y vector to 1
			for(int i = 0 ; i < n ; i++){
				x[i] = 1;
			}


			
			boolean flag = true;
			while(flag ){
				// set the flag to false assuming there is no convergence
				flag = false;
				// start computing the values of Y matrix
				
				for(int i = 0 ; i < n ; i++){
					double  x1= 0 ,x2 = 0;
					double temp = 0;
					for( int j = 0 ; j < i ; j++){
						x1+= x[j]*A[i][j];
					}
					for( int j = i+1 ; j <= n-1 ; j++){
						x2+= x[j]*A[i][j];
					}

					// compute new y value
					y[i] = (b[i] - x1 - x2)/A[i][i];
					// swap values
					}
					// swap x and y
				double [] temp = new double[n];
				temp = x;
				x=y;
				y = temp;
					//test convergence 
					for(int i = 0 ; i < n ; i++)
						if(Math.abs((2.0*(x[i] - y[i]))/(x[i] + y[i]))  >= (0.00000001))
					{
						flag = true;

					}


				
				// store the new value of y as pre-condition x
				
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		// print out the output

		if (n <= 100)
			for (int i = 0; i < n; ++ i)
				System.out.printf ("%d %g%n", i, y[i]);
		else
		{
			for (int i = 0; i <= 49; ++ i)
				System.out.printf ("%d %g%n", i, y[i]);
			for (int i = n - 50; i < n; ++ i)
				System.out.printf ("%d %g%n", i, y[i]);
		}
		long t2 = System.currentTimeMillis();

		
		System.out.printf ("%d msec%n", t2 - t1); 
	}
}

