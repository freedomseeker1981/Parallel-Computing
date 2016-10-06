import edu.rit.mp.BooleanBuf;
import edu.rit.mp.Buf;
import edu.rit.mp.DoubleBuf;

import edu.rit.mp.buf.BooleanItemBuf;

import edu.rit.pj.Comm;

import edu.rit.pj.reduction.BooleanOp;
/*
 * Name : Reda Elbahi

 * Course : Parallel computing I
 * Prof: Alan Kaminsky
 * Project: Programming Project I parallel version
 * Description: This program implements the jacobi method to solve linear algebra equations and test the convergnces
 * against the value of 10^-8 and stops the iteration at that level and prints the results .
 * The program uses the Cluster technique to solve the problem and find the final result in the process 0 and print it. 
 */
import edu.rit.util.Random;
import edu.rit.util.Range;


public class JacobiClu {

	static int n = 0;
	static double[][] A ;//={{8.0,2.0,3.0},{-1.0 ,6.0 ,4.0},{5.0 ,1.0, 7.0}} ;
	static double[] b ;//= {-4,5,9};

	static double[] x;//pre-condition
	static double[] y;// post-condition

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// slice variables
	static int ub;
	static int lb;
	static Range myslice;
	static Range[] slices;
	static int mylen;

	// convergence check variables
	static BooleanItemBuf outBuf ;
	static BooleanItemBuf inBuf;


	static DoubleBuf mOutbuf = null;
	static DoubleBuf[] mOutbufs = null;
	public static void main(String args[]) throws Exception{

		// Start timing.
		long t1 = System.currentTimeMillis();

		// Initialize world communicator.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();
		//	 long seed = Random seed
		if(args.length != 2)
		{
			System.err.println("worng input :   Please check your argument ");
			System.exit(1);
		}


		// slice the range 


		try{
			long seed = Long.parseLong(args[1]);
			n = Integer.parseInt(args[0]);


			slices = new Range(0,n-1) .subranges(size);
			myslice = slices[rank];
			lb = myslice.lb();
			ub = myslice.ub();
			mylen = myslice.length();
			// initialize buffers
			mOutbufs= new DoubleBuf[n];
			for(int i = 0 ; i < n ; i++)
				mOutbufs[i] = DoubleBuf.buffer();
			mOutbuf = DoubleBuf.buffer();

			// intialize the x and y value
			x = new double [n];
			y = new double [n];
			A = new double [mylen] [n];

			b = new double [mylen];

			//			mOutbuf = DoubleBuf.sliceBuffer(y,myslice);
			Random prngTh = Random.getInstance (seed);

//			System.out.println(" lower bound " + lb + " upper bound " + ub);
			prngTh.setSeed(seed);
			prngTh.skip(lb*(n + 1));
			for (int i = 0; i < mylen; ++ i)
			{
				for (int j = 0; j < n; ++ j)
					A[i][j] = prngTh.nextDouble()*9.0 + 1.0;
				A[i][i+lb] += 10.0*n;// jump by lower bound to prevent overlapping
				b[i] = prngTh.nextDouble()*9.0 + 1.0;
			}

			for(int i = 0 ; i < n ; i++)
				x[i] = 1;
			inBuf = new BooleanItemBuf(false);

			while(inBuf.item == false){
				inBuf.item = true;

				// Calculate the new value
				for(int i = 0; i < mylen ; i++){
					double  x1= 0.0 ,x2 = 0.0;

					for( int j = 0; j < i+lb ; j++){
						x1+= x[j]*A[i][j];
					}
					for( int j = i+1+lb ; j < n ; j++){
						x2+= x[j]*A[i][j];
					}

					// compute new y value
					y[i+ lb] = (b[i] - x1 - x2)/A[i][i+lb];// jump to i + lower bound


				}// end of inner for loop

				mOutbuf = DoubleBuf.sliceBuffer(y,myslice);
				world.allGather(mOutbuf,mOutbufs);
				world.barrier();
				// swap the x and y
				double[] temp = new double[n];
				temp = y;
				y = x ;
				x = temp;
				//System.out.println(" swapping");

				boolean tempflag = true;


				for(int i =lb ; i <= ub ; i++)
					if(Math.abs((2.0*(x[i] - y[i]))/(x[i] + y[i]))  >= (0.00000001))
					{
						// temp flag for each thread will affect the while loop flag
						tempflag = false;

					}


				// check for all processes the result of the convergence
				inBuf.item = tempflag;
				world.allReduce (inBuf, BooleanOp.AND);

			}// end of while loop
		}catch(Exception e){
			e.printStackTrace();
			System.exit(1);
		}


		// gather all portions into process 0 
		//		mOutbufs = DoubleBuf.sliceBuffers( y , slices);


		// print final result for only process 
		if(rank == 0){
			if (n <= 100)
				for (int i = 0; i < n; ++ i)
					System.out.printf ("%d %g%n", i, x[i]);
			else
			{
				for (int i = 0; i <= 49; ++ i)
					System.out.printf ("%d %g%n", i,y[i] );
				for (int i = n-50 ; i < n; ++ i)
					System.out.printf ("%d %g%n", i,y[i]);
			}
			long t2 = System.currentTimeMillis();
			System.out.printf ("%d msec%n", t2-t1); 

		}
	}
}


