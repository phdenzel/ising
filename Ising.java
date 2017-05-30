/*
   @author: dephil

   Simple 2D transverse field Ising model solver using MCMC

   Compilation (command line):
     > javac Ising.java
     > java Ising <Args>

   Args:
     mode  -  (String)  decides in which mode the program runs ('gx' or 'plot')
     L     -  (Integer) lattice dimension (e.g. 100)
     T     -  (Double)  temperature (value between 1 to 4 for instance; critical temperature theoretically at ~2.26918; not needed for 'plot')

   Usage (command line):
     > java Ising gx 100 2.2
     or
     > java Ising plot 100

   Note: Makes use of the jfreechart-1.0.17 package.
         Download: http://www.jfree.org/jfreechart/download.html
         For LINUX/MacOS add the path of the necessary .jar files in jfreechart-1.0.17/lib/ to the CLASSPATH in your .bashrc/.profile file; for example:
         > export CLASSPATH=$CLASSPATH:/home/dephil/local/jfreechart-1.0.17/lib/jfreechart-1.0.17.jar:/home/dephil/local/jfreechart-1.0.17/lib/jcommon-1.0.21.jar
         For WINDOWS add the following two Strings to the end of the user variable CLASSPATH:
         > "C:\jfreechart-1.0.17\lib\jfreechart-1.0.17.jar"
         > "C:\jfreechart-1.0.17\lib\jcommon-1.0.21.jar"
 */
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


public class Ising extends JPanel {

// Model content
int L;                // linear dimension of lattice
int N;                // number of spins
boolean[][] spin;     // 2d lattice of spin sites (true=up, false=down)

double J = 1.;        // coupling constant (ferromagnetic for J>0)
double B = 0.;        // transverse field

double T;             // temperature (kB=1)
double M;             // total magnetic moment
double H;             // Hamiltonian aka total energy

double[][] bF;        // Boltzmann factors

int MCsteps;          // total Monte-Carlo steps
int flipsteps;        // accepted Monte-Carlo steps

double Msum;          // sum of magnetic moment
double M2sum;         // quadratic sum of magnetic moment
double Hsum;          // sum of total energy
double H2sum;         // quadratic sum of total energy

double MAvg;          // average magnetic moment
double HAvg;          // average energy
double Cv;            // average specific heat
double chi;           // average magnetic susceptibility
double flipAvg;       // spin flip average

// JPanel content
final static int XWIN = 600;     // graphics window width
final static int YWIN = 600;     // graphics window height
final static int XLOC = 100;     // x distance from upper left corner on screen
final static int YLOC = 100;     // y distance from upper left corner on screen
int PIXEL_SIZE;                  // pixel size is determined by the lattice and window size
private JFreeChart chart;

private boolean plot_mode;       // true if plot should be shown
/* Ising model ************************************************************** */

public Ising(int L, double T, double prob) { // prob = probability of spin up at each site
        this.L = L;
        this.N = L*L;
        this.T = T;
        this.spin = new boolean[L][L];
        // initialize random spins and calculate initial magnetization / energy
        int m = 0;
        double e = 0.;
        int sm = 0;
        for (int i=0; i<L; i++) {
                // neighbor indices (periodic boundary, hence +L)
                int east = i+1;
                int west = i-1+L;
                for (int j=0; j<L; j++) {
                        // random spin assignment according to prob
                        spin[i][j] = (Math.random() < prob);
                        // magnetization
                        if (spin[i][j])
                                m++;
                        else
                                m--;
                        // staggered magnetization
                        if ((i+j) % 2 == 0) {
                                if (spin[i][j]) sm++;
                                else sm--;
                        }
                        else {
                                if (spin[i][j]) sm--;
                                else sm++;
                        }
                        // energy
                        if (spin[i][j] == spin[(east)%L][j]) e++;
                        else e--;
                        if (spin[i][j] == spin[i][(j+1)%L]) e++;
                        else e--;
                        if (spin[i][j] == spin[(west)%L][j]) e++;
                        else e--;
                        if (spin[i][j] == spin[i][(j-1+L)%L]) e++;
                        else e--;
                }
        }
        // initial magnetic moment and energy
        this.M = (J>=0) ? (double)m : (double)sm;
        this.H = -.5*J*e + B*(double)m; // correct for double counting
        // initialize Monte-Carlo variables
        this.MCsteps = 0;
        this.flipsteps = 0;
        // initialize rest of the variables
        initialize_sums();
        initialize_avgs();
        // precalculate Boltzmann factors
        calcBoltzmann();
}


void calcBoltzmann() {
        // create precomputed table of Boltzmann factors
        if (bF == null)
                bF = new double[17][3];

        for (int i=-8; i<=8; i+=4) {
                bF[i+8][0] = Math.exp(-(i*J+2*B)/T);
                bF[i+8][2] = Math.exp(-(i*J-2*B)/T);
        }
}

public double magnetization() {
        int m = 0;
        int sm = 0;
        for (int i=0; i<L; i++) {
                for (int j=0; j<L; j++) {
                        if (spin[i][j])
                                m++;
                        else
                                m--;
                        if ((i+j) % 2 == 0) {
                                if (spin[i][j]) sm++;
                                else sm--;
                        }
                        else {
                                if (spin[i][j]) sm--;
                                else sm++;
                        }
                }
        }
        // return total magnetic moment depending on coupling constant
        return (J>=0) ? (double)m : (double)sm;
}

public int site_interaction(int n, int m) {
        int snsm = 0;
        if (spin[n][m] == spin[(n+1)%L][m]) snsm++;
        else snsm--;
        if (spin[n][m] == spin[n][(m+1)%L]) snsm++;
        else snsm--;
        if (spin[n][m] == spin[(n-1+L)%L][m]) snsm++;
        else snsm--;
        if (spin[n][m] == spin[n][(m-1+L)%L]) snsm++;
        else snsm--;
        return snsm;
}

public double energy() {
        double m = 0.;
        double e = 0.;
        for (int i=0; i<L; i++) {
                for (int j=0; j<L; j++) {
                        e += site_interaction(i, j);
                        if (spin[i][j])
                                m++;
                        else
                                m--;
                }
        }
        return -.5*J*e + B*m;   // correct for double counting by dividing by 2
}

public void update_M() {
        M = magnetization();
}

public void update_H() {
        H = energy();
}

public void metropolis() {
        int dsm = 0;
        for (int si=0; si<N; si++) {
                // pick random spin site
                int x = (int)(L*Math.random());
                int y = (int)(L*Math.random());
                int spinxy = (spin[x][y]) ? 1 : -1;
                int delEpJ = 2*site_interaction(x, y);
                double delE = delEpJ*J;
                if ((Math.random() < bF[delEpJ+8][spinxy+1]) || (delE <= 0)) {
                        //if ((Math.random() < Math.exp(-delE/T)) || (delE <= 0)) {
                        // flip spin
                        spin[x][y] = !spin[x][y];
                        ++flipsteps;
                        // update total magnetic moment / energy
                        if ((x+y) % 2 == 0) {
                                if (spin[x][y]) dsm++;
                                else dsm--;
                        }
                        else {
                                if (spin[x][y]) dsm--;
                                else dsm++;
                        }
                        M += (J>=0) ? (double)(2*spinxy) : (double)dsm;
                        H += delE;
                }
        }
        ++MCsteps;
        // update total magnetic moment / energy
        // update_M();
        // update_H();
        // update summed quantities
        update_sums();
        // update averaged quantities
        update_avgs();
}

public void initialize_sums() {
        this.Msum = 0;
        this.M2sum = 0;
        this.Hsum = 0;
        this.H2sum = 0;
}

public void initialize_avgs() {
        this.flipAvg = 0;
        this.HAvg = 0;
        this.MAvg = 0;
        this.Cv = 0;
        this.chi = 0;
}

public void update_sums() {
        Msum += M;
        M2sum += M*M;
        Hsum += H;
        H2sum += H*H;
}

public void update_avgs() {
        double invN = 1./N;
        if (MCsteps > 0)
                invN /= MCsteps;
        flipAvg = flipsteps*invN;
        MAvg = Msum*invN;
        HAvg = Hsum*invN;
        Cv = (H2sum*invN - N*HAvg*HAvg) / (T*T);
        chi = (M2sum*invN - N*MAvg*MAvg) / T;
}

public void set_temperature(double T) {
        this.T = T;
}

/* JPanel ******************************************************************* */

public void pixel_size() {
        this.PIXEL_SIZE = XWIN/L; // if integer division is not exact the window will be too big
}

public void draw_gx(Graphics g) {
        Rectangle bounds = getBounds();
        g.clearRect(bounds.x, bounds.y, bounds.width, bounds.height);

        /* Draw the grid. */
        for (int i=0; i<L; i++) {
                for (int j=0; j<L; j++) {
                        if (spin[i][j])
                                g.setColor(Color.WHITE);
                        else
                                g.setColor(Color.GRAY);
                        int x = PIXEL_SIZE * i;
                        int y = PIXEL_SIZE * j;
                        g.fillRect(x, y, PIXEL_SIZE, PIXEL_SIZE);
                }
        }
}

public void draw_plot(Graphics g) {
        if (chart != null) {
                Rectangle b = g.getClipBounds();
                chart.draw((Graphics2D)g, new Rectangle2D.Double(b.x, b.y, b.width, b.height));
        }
}

public void set_plot_mode(boolean plot) {
        this.plot_mode = plot;
}

public void paint(Graphics g) {
        super.paint(g);
        if (plot_mode) {
                draw_plot(g);
        }
        else {
                draw_gx(g);
        }
}

public void run_gx() {
        // assumes that an Ising object was already initialized
        pixel_size();
        set_plot_mode(false);
        while (true) { // runs until window is closed
                metropolis();
                repaint();
                try {
                        Thread.sleep(50);
                } catch (InterruptedException e) {
                        e.printStackTrace();
                }
        }
}

public void run_plot() {
        // assumes that an Ising object was already initialized
        set_plot_mode(true);

        // allocate data arrays
        int d = 100; // array size
        double T_o = 0.1; // start measurement
        double T_d = 4.; // end measurement
        double dT = (T_d-T_o) / d;
        double[] temp = new double[d];
        double[] avgm = new double[d];
        double[] avge = new double[d];
        double[] avgcv = new double[d];
        double[] avgchi = new double[d];
        double[] fliprate = new double[d];

        // initialize temperature data
        temp[0] = T_o;
        for (int i=1; i<d; i++) {
                temp[i] = temp[i-1]+dT;
        }

        // produce data
        for (int i=0; i<d; i++) {
                set_temperature(temp[i]); // set temperature
                calcBoltzmann(); // redo Boltzmann table
                //for (int j=0; j<1; j++)
                metropolis();
                avgm[i] = MAvg;
                avge[i] = HAvg;
                avgcv[i] = Cv;
                avgchi[i] = chi;
                fliprate[i] = flipAvg;
        }

        // for JFreeChart
        XYSeries series = new XYSeries("");
        // XYSeries series = new XYSeries("total energy");
        for (int i=0; i<d; i++) {
                series.add(temp[i], avge[i]);
        }
        XYDataset dataset = new XYSeriesCollection(series);
        chart = ChartFactory.createXYLineChart("Ising 2D model", "temperature", "average total energy", dataset);

        // plot
        repaint();
}

/* Main method ************************************************************** */

public static void main(String[] args) {
        // parse command line args
        String mode = ""; // mode arg
        if(args.length==0) {
                System.out.println("No arguments given...");
                System.out.println("Will run in 'gx' mode with a 100x100 lattice at temperature 2.0");
                System.out.println("Rerun with, e.g.");
                System.out.println("> java Ising gx 200 1.5");
                System.out.println("> java Ising plot 100");
                mode = "gx";
        }
        else {
                mode = args[0];
        }
        boolean gx = mode.equals("gx");
        boolean plot = mode.equals("plot");

        int n = 100; // lattice dimension arg
        if(args.length>1) {
                n = Integer.parseInt(args[1]);
        }

        double t = 0.1; // temperature arg
        if (gx && args.length>2) {
                t = Double.parseDouble(args[2]);
        }

        // create an Ising object
        if (gx || plot) {
                Ising ising = new Ising(n, t, 0.5);  // default prob=0.5
                JFrame top = new JFrame("2D Ising model");
                top.setBounds(XLOC, YLOC, XWIN, YWIN);
                top.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                top.getContentPane().add(ising);
                top.setVisible(true);
                // graphics
                if (gx) {
                        ising.run_gx();
                }
                // plot
                else {
                        ising.run_plot();
                }
        }
        else {
                System.out.println("Problem detected! Exiting...");
                System.exit(0);
        }
}

} /* END ISING ********************************************************* */
